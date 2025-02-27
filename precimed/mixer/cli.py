#!/usr/bin/env python

import argparse
import collections
import json
import numpy as np
import os
import pandas as pd
import random
import scipy.optimize
import scipy.stats
import scipy.io as sio
import time
import sys
import intervaltree

from datetime import datetime
from pathlib import Path
from scipy.sparse import csr_matrix, coo_matrix

from .libbgmg import LibBgmg, LibBgmgAnnot
from .utils import AnnotUnivariateParams
from .utils import BivariateParams
from .utils import _params_to_dict
from .utils import _dict_to_params
from .utils import UnivariateParametrization_natural_axis		                     # diffevo, neldermead, uncertainty
from .utils import UnivariateParametrization_constPI_constSIG2BETA                   # inflation
from .utils import UnivariateParametrization_constPI		                         # infinitesimal
from .utils import BivariateParametrization_constUNIVARIATE_constRHOBETA_constPI     # inflation
from .utils import BivariateParametrization_constUNIVARIATE_natural_axis             # diffevo, neldermead
from .utils import BivariateParametrization_constUNIVARIATE_constRG_constRHOZERO     # brute1, brent1
from .utils import BivariateParametrization_constUNIVARIATE_infinitesimal            # infinitesimal
from .utils import AnnotUnivariateParametrization

from .utils import _log_exp_converter

from .utils import _calculate_univariate_uncertainty
from .utils import _calculate_univariate_uncertainty_funcs
from .utils import _calculate_bivariate_uncertainty_funcs

from .utils_obsolete import UnivariateParams_obsolete

from .utils import calc_qq_plot
from .utils import calc_qq_plot_plsa
from .utils import calc_power_curve
from .utils import calc_bivariate_qq
from .utils import calc_bivariate_pdf

from .utils import NumpyEncoder

from .utils import _auxoption_tagpdf
from .utils import _auxoption_none

MASTHEAD = "***********************************************************************\n"
MASTHEAD += "* (c) 2016-2024 MiXeR software - Univariate and Bivariate Causal Mixture for GWAS\n"
MASTHEAD += "* Norwegian Centre for Mental Disorders Research / University of Oslo\n"
MASTHEAD += "* Center for Multimodal Imaging and Genetics / UCSD\n"
MASTHEAD += "* GNU General Public License v3\n"
MASTHEAD += "***********************************************************************\n"

_cost_calculator_sampling = 0
_cost_calculator_gaussian = 1
_cost_calculator_convolve = 2
_cost_calculator_smplfast = 3

_cost_calculator = {
    'gaussian': _cost_calculator_gaussian, 
    'sampling': _cost_calculator_sampling,
    'convolve': _cost_calculator_convolve,
    'smplfast': _cost_calculator_smplfast,
}

GO_EXTEND_BP_DEFAULT = 10000

CHR_OFFSET = np.int64(1e10)

def enhance_optimize_result(r, cost_n, cost_df=None, cost_fast=None):
    # optimize_result - an instance of scipy.optimize.OptimizeResult
    # const_n - number of genetic variants effectively contributing to the cost (sum of weights)
    r['cost_n'] = cost_n
    r['cost_df'] = len(r.x) if (cost_df is None) else cost_df
    r['cost'] = r.fun
    r['BIC'] = np.log(r['cost_n']) * r['cost_df'] + 2 * r['cost']
    r['AIC'] =                   2 * r['cost_df'] + 2 * r['cost']
    if cost_fast is not None: r['cost_fast'] = cost_fast

def check_input_file(args, argname, chri=None):
    arg_dict = vars(args)
    if (argname in arg_dict) and arg_dict[argname]:
        argval = arg_dict[argname]
        if chri: argval=argval.replace('@', str(chri))
        if not os.path.isfile(argval): raise ValueError("Input file --{a} does not exist: {f}".format(a=argname, f=argval))

def fix_chr2use_arg(args, argname="chr2use"):
    arg_dict = vars(args)
    if argname not in arg_dict : return
    if not arg_dict[argname]: arg_dict[argname] = []; return
    chr2use_arg = arg_dict[argname]
    chr2use = []
    if chr2use_arg is not None:
        for a in chr2use_arg.split(","):
            if "-" in a:
                start, end = [int(x) for x in a.split("-")]
                chr2use += [str(x) for x in range(start, end+1)]
            else:
                chr2use.append(a.strip())
        if np.any([not x.isdigit() for x in chr2use]): raise ValueError('Chromosome labels must be integer')
    arg_dict[argname] = [int(x) for x in chr2use]

def make_ranges(args_exclude_ranges):
    # Interpret --exclude-ranges input
    ChromosomeRange = collections.namedtuple('ChromosomeRange', ['chr', 'from_bp', 'to_bp'])
    exclude_ranges = []
    if args_exclude_ranges is not None:
        for exclude_range in args_exclude_ranges:
            try:
                exclude_range=exclude_range.strip().upper().replace('CHR','')
                if exclude_range == '[]': continue
                if exclude_range == 'MHC': exclude_range='6:25-35MB'
                if exclude_range == 'APOE': exclude_range = '19:45-46MB'
                if exclude_range.endswith('KB'): factor=1e3
                elif exclude_range.endswith('MB'): factor=1e6
                else: factor=1
                exclude_range=exclude_range.replace('KB', '').replace('MB', '')
                chr_from_to = [int(x) for x in exclude_range.replace(':', ' ').replace('-', ' ').split()[:3]]
                chr_from_to[1] *= factor
                chr_from_to[2] *= factor
                range = ChromosomeRange._make(chr_from_to)
            except Exception as e:
                raise(ValueError('Unable to interpret exclude range "{}", chr:from-to format is expected.'.format(exclude_range)))
            exclude_ranges.append(range)
    return exclude_ranges

def fix_and_validate_args(args):
    # here is what you need to know about how --seed  is handled:
    # 1. If it optional. When not specified, we will use np.rando.randint to generate a random seed (the the line below this comment)
    # 2. There are several places in this code where we generate random numbers using numpy library. Therefore, we use args.seed to initialize np.random.seed.
    # 3. We pass args.seed to MiXeR using set_options("seed", args.seed), see convert_args_to_libbgmg_options() function.
    #    This will ensure that random prunning uses consistent seed.
    # 4. In ADAM algorithm, each cost function evaluation starts with it's own seed (to respect the stochastic nature of the algorithm)
    #    This is done in line "libbgmg.set_option('seed', np.random.randint(np.iinfo(np.int32).max))".
    # 5. --diffevo-fast-repeats also rely on seed to ensure different path of the differential evolution algorihtm.
    if 'seed' in args:
        if args.seed is None: args.seed = np.random.randint(np.iinfo(np.int32).max)
        np.random.seed(args.seed)

    fix_chr2use_arg(args, argname="chr2use")
    if 'trait1_file' in args:
        if '@' in args.trait1_file:
            for chri in args.chr2use:
                    check_input_file(args, 'trait1_file', chri)
        else:
            check_input_file(args, 'trait1_file')
    check_input_file(args, 'trait2_file')
    check_input_file(args, 'trait1_params_file')
    check_input_file(args, 'trait2_params_file')
    check_input_file(args, 'load-params-file')

    if ('analysis' in args) and (args.analysis == 'fit2'):
        if not args.trait1_params_file: raise ValueError('--trait1-params-file is required ')
        if not args.trait2_params_file: raise ValueError('--trait2-params-file is required ')

    arg_dict = vars(args)
    for chri in arg_dict["chr2use"]:
        check_input_file(args, 'bim-file', chri)
        check_input_file(args, 'ld-file', chri)
        check_input_file(args, 'annot-file', chri)

def convert_args_to_libbgmg_options(args):
    libbgmg_options = {
        'r2min': args.r2min if ('r2min' in args) else None,
        'kmax': args.kmax[0] if ('kmax' in args) else None, 
        'kmax_pdf': args.kmax_pdf[0] if ('kmax_pdf' in args) else None,
        'threads': args.threads[0] if ('threads' in args) else None,
        'seed': args.seed if ('seed' in args) else None,
        'cubature_rel_error': args.cubature_rel_error if ('cubature_rel_error' in args) else None,
        'cubature_max_evals': args.cubature_max_evals if ('cubature_max_evals' in args) else None,
        'z1max': args.z1max if ('z1max' in args) else None,
        'z2max': args.z2max if ('z2max' in args) else None, 
    }
    return [(k, v) for k, v in libbgmg_options.items() if v is not None ]

def set_weights(args, libbgmg, randprune_n):
    if 'weights_file' in args:
        if args.weights_file.lower() == 'auto':
            if ('load_baseline_params_file' in args) and args.load_baseline_params_file:
                args.weights_file = (args.load_baseline_params_file + '\n').replace('.json\n', '.weights', 1)
            elif ('load_params_file' in args) and args.load_params_file:
                args.weights_file = (args.load_params_file + '\n').replace('.json\n', '.weights', 1)
            else:
                args.weights_file = 'none'

        if args.weights_file.lower() != 'none':
            args.save_weights = False
            check_input_file(args, 'weights_file')

            df_weights = pd.read_csv(args.weights_file, sep='\t')
            df_weights = df_weights[df_weights['chrnumvec'].isin(set(libbgmg.chrnumvec))]
            if len(df_weights) != libbgmg.num_snp:
                raise ValueError(f'--weights-file {args.weights_file} defines weights for {len(df_weights)} snps, while {libbgmg.num_snp} snps were expected')
            libbgmg.log_message(f'apply --weights-file {args.weights_file} (and skip all --hardprune and --randprune settings)')
            libbgmg.weights = df_weights.loc[libbgmg.defvec, 'weights'].values.flatten()
            return

    hardprune_subset = args.hardprune_subset if ('hardprune_subset' in args) else 0
    hardprune_r2 = args.hardprune_r2 if ('hardprune_r2' in args) else 0
    hardprune_maf = args.hardprune_maf if ('hardprune_maf' in args) else 0

    if (hardprune_subset != 0) or (hardprune_r2 != 0) or (hardprune_maf != 0):
        libbgmg.set_weights_hardprune(hardprune_subset, hardprune_r2, hardprune_maf, not args.disable_inverse_ld_score_weights)

    if (randprune_n is not None) and (randprune_n != 0):
        libbgmg.set_weights_randprune(randprune_n, args.randprune_r2, 0, not args.disable_inverse_ld_score_weights, "", "")

# https://stackoverflow.com/questions/27433316/how-to-get-argparse-to-read-arguments-from-a-file-with-an-option-rather-than-pre
class LoadFromFile (argparse.Action):
    def __call__ (self, parser, namespace, values, option_string=None):
        with values as f:
            contents = '\n'.join([x for x in f.readlines() if (not x.strip().startswith('#'))])

        data = parser.parse_args(contents.split(), namespace=namespace)
        for k, v in vars(data).items():
            if v and k != option_string.lstrip('-'):
                setattr(namespace, k, v)

def parser_add_arguments(parser, func, analysis_type):
    if analysis_type in ['split_sumstats']:
        parser.add_argument("--trait1-file", type=str, default="", help="GWAS summary statistics for the first trait (required argument);")

    if analysis_type in ['fit1', 'test1', 'fit2', 'test2', 'plsa', 'perf', 'fixed_effects', 'snps']:
        parser.add_argument("--bim-file", type=str, default=None, help="plink bim file (required argument); "
            "defines the reference set of SNPs used for the analysis. "
            "Marker names must not have duplicated entries. "
            "May contain symbol '@', which will be replaced by an actual chromosome label. ")
        parser.add_argument("--ld-file", type=str, default=None, help="file with linkage disequilibrium information, "
            "generated via 'mixer.py ld' command (required argument); "
            "may contain symbol '@', similarly to --bim-file argument. ")

    if True:  # all analysis types
        parser.add_argument("--chr2use", type=str, default="1-22", help="chromosome labels to use (default: %(default)s); "
        "chromosome must be labeled by integer, i.e. X and Y are not acceptable; example of valid arguments: '1,2,3' or '1-4,12,16-20'")

    if analysis_type in ['fit1', 'test1', 'fit2', 'test2', 'plsa', 'perf', 'fixed_effects']:
        parser.add_argument('--extract', type=str, default="", help="(optional) File with variants to include in the analysis. By default, all variants are included. This applies to GWAS tag SNPs, however LD is still computed towards the entire reference provided by --bim-file. This applies before --exclude, so a variant listed in --exclude won't be re-introduced if it's also present in --extract list.")
        parser.add_argument('--exclude', type=str, default="", help="(optional) File with variants to exclude from the analysis.")
        parser.add_argument('--exclude-ranges', type=str, default=['MHC'], nargs='+',
            help='(default: %(default)s) exclude SNPs in ranges of base pair position; '
            'the syntax is chr:from-to, for example 6:25000000-35000000; '
            'multiple regions can be excluded; '
            '"chr" prefix prior to chromosome label, as well as KB and MB suffices are allowed, e.g. chr6:25-35MB is a valid exclusion range. '
            'Some special case regions are also supported, for example "--exclude-ranges MHC APOE".'
            'To overwrite the default, pass "--exclude-ranges []".')

        parser.add_argument('--use-complete-tag-indices', default=False, action="store_true", help=argparse.SUPPRESS) #"advanced option (expert use only); MiXeR software has a notion of 'tag SNPs', which are the subset of the SNPs in --bim-file with well defined GWAS z-score(s); this saves memory used by LD r2 storage, which became a rectangular matrix (LD r2 between two SNPs from the --bim-file is stored only if one of the SNPs is a tag SNP). Setting this flag enforces all --bim-file SNPs to be treated as tag SNPs, some will have undefined GWAS z-score(s), and the complete LD matrix will be save in memory.")
        parser.add_argument('--allow-ambiguous-snps', default=False, action="store_true", help="advanced option (expert use only); a flag allowing to include A/T and C/G SNPs in fit procedure.")
        parser.add_argument("--savelib-file", type=str, default=None, help="optional path to save the state of the libbgmg object (for faster initialization with --loadlib-file). It is highly recommended to use --use-complete-tag-indices together with --savelib-file.")
        parser.add_argument("--loadlib-file", type=str, default=None, help="optional path to load libbgmg object; "
            "there are several limitations if you use this option: "
            "--bim-file argument must be the same as what was used with --savelib-file; "
            "same restriction applies to --chr2use argument for fit1,test1,fit2,test2 analyses, however for plsa analysis --chr2use can be different as each chromosome has its own .bin file; "
            "--ld-file and --use-complete-tag-indices arguments will be ignored.")

    if analysis_type in ['fit1', 'test1', 'fit2', 'test2', 'plsa', 'perf']:
        num_traits = None
        if analysis_type in ['fit1', 'test1', 'plsa']: num_traits = 1
        if analysis_type in ['fit2', 'test2', 'perf']: num_traits = 2
        is_fit = analysis_type in ['plsa', 'fit1', 'fit2']

        parser.add_argument("--trait1-file", type=str, default="", help="GWAS summary statistics for the first trait (required argument); for 'plsa' analysis it is recommended to split GWAS summary statistics per chromosome; if this is done then --trait1-file should contain symbol '@', which will be replaced by an actual chromosome label. ")
        if num_traits==2: parser.add_argument("--trait2-file", type=str, default="", help="GWAS summary statistics for the second trait (required argument)")
        parser.add_argument('--z1max', type=float, default=None, help="right-censoring threshold for the first trait (default: %(default)s); recommended setting: '--z1max 9.336' (equivalent to p-value 1e-20)")
        if num_traits==2: parser.add_argument('--z2max', type=float, default=None, help="right-censoring threshold for the second trait (default: %(default)s); recommended setting: '--z2max 9.336'")
        parser.add_argument('--hardprune-subset', type=float, default=0, help="fraction of SNPs to randomly select (default: %(default)s); passing 0 as argument disables the operation; values below 1 are interpreted as fraction, values above 1 are interpreted as number of SNPs to select. Both options consider all SNPs in the reference (not just tag SNPs)")
        parser.add_argument('--hardprune-maf', type=float, default=0.05, help="threshold for minor allele frequency (default: %(default)s); passing 0 as argument disables the operation")
        parser.add_argument('--hardprune-r2', type=float, default=1.0, help="LD r2 threshold for prunning SNPs (default: %(default)s); passing 0 as argument disables the operation")
        parser.add_argument('--randprune-n', type=int, default=0, help="number of random pruning iterations (default: %(default)s)")
        parser.add_argument('--randprune-r2', type=float, default=0.1, help="threshold for random pruning (default: %(default)s)")
        parser.add_argument('--disable-inverse-ld-score-weights', default=False, action="store_true", help="a flag allowing to disable weighting by inverse ld score")

        parser.add_argument('--weights-file', type=str, default=('auto' if (analysis_type in ['fit1', 'plsa']) else 'none'), help="file to load weights (default: %(default)s); 'auto' takes weights file that is associated with --load-params argument (if it is provided); 'none' disables 'auto' feature; otherwise this argument must point to a <out>.weights produced by a previous mixer.py analysis")
        parser.add_argument('--save-weights', default=(analysis_type == 'plsa'), action="store_true", help="a flag allowing to save weights set via--hardprune and --randprune options; does not have effect is weights are loaded from --weights-file")

        parser.add_argument('--seed', type=np.uint32, default=None, help="random seed (default: %(default)s).")
        parser.add_argument('--r2min', type=float, default=0.0, help=argparse.SUPPRESS) # help="r2 values below this threshold will contribute via infinitesimal model")
        parser.add_argument('--threads', type=int, default=[None], nargs='+', help="specify how many concurrent threads to use for computations; (default: total number of CPU cores available)")

        if analysis_type == 'perf':
            parser.add_argument('--kmax', type=int, default=[20000, 2000, 200], nargs='+', help="number of sampling iterations")
        else:
            parser.add_argument('--kmax', type=int, default=[20000], nargs='+', help="number of sampling iterations for log-likelihod and posterior delta (default: 20000)")
            parser.add_argument('--kmax-pdf', type=int, default=[100], nargs='+', help=argparse.SUPPRESS) # help="number of sampling iterations for qq-plot and power-curves (default: 100)")

    if analysis_type == 'fixed_effects':
        parser.add_argument("--causal-variants", type=str, default=None, help="File with causal effects (similar to --causal-variants in https://github.com/precimed/simu)")

    if analysis_type == 'plsa':
        # --go-file and --go-file-test must contain columns 'GO', 'GENE', 'CHR', 'FROM', 'TO'
        parser.add_argument("--annot-file", type=str, default=None, help="(optional) path to binary annotations in LD score regression format, i.e. <path>/baseline.@.annot.gz for fitting enrichment model model. This must include the first column with all ones ('base' annotation category covering the entire genome).")
        parser.add_argument("--annot-file-test", type=str, default=None, help="(optional) path to binary annotations in LD score regression format, i.e. <path>/baseline.@.annot.gz for evaluating enrichment.  If provided, enrichment scores computed for --annot-file-test will be saved to <out>.annot_test_enrich.csv.")
        parser.add_argument("--go-file", type=str, default=None, help="(optional) path to GO antology file for fitting enrichment model model. The format is described in the documentation. 'base' category that covers entire genome will be added automatically.")
        parser.add_argument("--go-file-test", type=str, default=None, help="(optional) path to an additional  GO antology file for evaluating enrichment, same convention as --go-file. If provided, enrichment scores computed for --go-test-file will be saved to <out>.go_test_enrich.csv")
        parser.add_argument("--calc-loglike-diff-go-test", default=None, choices=['fast', 'full'], help="algorithm for computation of the loglike_diff column in <out>.go_test_enrich.csv")
        parser.add_argument("--go-all-genes-label", type=str, default='coding_genes', help="reference gene-set to calibrate fold enrichment, e.g. allowing to compute enrichment w.r.t. the set of all coding genes (default: %(default)s)")
        parser.add_argument("--go-extend-bp", type=int, default=GO_EXTEND_BP_DEFAULT, help="extends each gene by this many base pairs, defining a symmetric window up and downstream (default: %(default)s)")
        parser.add_argument('--save-ldsc-reference', default=False, action="store_true", help=argparse.SUPPRESS) #"advanced option (expert use only); saves --annot-file and --go-file as LDSC reference")
        parser.add_argument('--enable-faster-gradient-calculation', default=False, action="store_true", help=argparse.SUPPRESS) #"advanced option (expert use only); saves --annot-file and --go-file as LDSC reference")
        
        parser.add_argument('--sig2-zeroL', default=0.0, type=float, help="(optional) initial value or constraint for 'sig2_zeroL' parameter (default: %(default)s); make sure to change initial value of this parameter to 1e-8 if '--fit sig2-zeroL'")
        parser.add_argument('--sig2-zeroA', default=1.0, type=float, help="(optional) initial value or constraint for 'sig2_zeroA' parameter (default: %(default)s);")
        parser.add_argument('--s-value', default=-0.25, type=float, help="(optional) initial value or constraint for the 's' parameter (default: %(default)s);")
        parser.add_argument('--l-value', default=0, type=float, help="(optional) initial value or constraint for the 'l' parameter (default: %(default)s);")
        parser.add_argument('--pi-value', default=1.0, type=float, help="(optional) initial value or constraint for the 'pi' parameter (default: %(default)s);")
        parser.add_argument('--h2-init-calibration', default=None, type=float, help="(optional) initial value for heritability calibration")
        parser.add_argument('--constrain-base-gene-category', default=False, action="store_true", help="a flag indicating whether to constrain sig2_gene for the baseline category during fit")
        parser.add_argument('--nullify-go-all-genes-sig2-gene', default=False, action="store_true", help="a flag indicating whether to set sig2_gene=0 for --go-all-genes-label")

        parser.add_argument('--fit', type=str, nargs='+', default=[], help="parameters to fit from the data",
            choices=['sig2-zeroL', 'sig2-zeroA', 's', 'l', 'pi', 'annot', 'gene'])

        # --fit annot gene s l pi sig2-zeroA sig2-zeroL 
        parser.add_argument('--annot-p', default=1, type=float, help="power factor for sigma2 aggregation in overlapping annotations (default: %(default)s)")
        parser.add_argument('--gene-p', default=1, type=float, help="power factor for sigma2 aggregation in overlapping gene-sets (default: %(default)s)")

        parser.add_argument('--adam-epoch', default=[10, 10, 10, 10, 10, 10, 10, 10, 10, 10], nargs='+', type=int, help="number of iterations in ADAM procedure (default: %(default)s)")
        parser.add_argument('--adam-beta1', default=0.9, type=float, help="beta_1 parameter in ADAM procedure (default: %(default)s)")
        parser.add_argument('--adam-beta2', default=0.99, type=float, help="beta_2 parameter in ADAM procedure (default: %(default)s)")
        parser.add_argument('--adam-eps', default=1e-8, type=float, help="epsilon parameter in ADAM procedure (default: %(default)s)")
        parser.add_argument('--adam-step', default=[0.064, 0.032, 0.016, 0.008, 0.004, 0.002, 0.001, 0.0005, 0.00025, 0.0001], nargs='+', type=float, help="step parameter in ADAM procedure (default: %(default)s)") # (two arguments trigger cycling learning rate")
        parser.add_argument('--adam-disable', default=False, action="store_true", help="a flag allowing to disable optimization; typical usecase would be in conjunction with these flags: '--adam-disable --load-params-file <out-of-a-previous-run>.json --make-snps-file --allow-ambiguous-snps'")
        parser.add_argument('--adam-separate-by-chromosome', default=False, action="store_true", help="a flag allowing to run optimization separately on each chromosome; this is only recommended with '--fit gene --constrain-base-gene-category --adam-beta2 0.9' configuration")
        
        parser.add_argument('--load-params-file', type=str, default=None, help="(optional) params of the fitted model to load; expected to be from a 'mixer.py plsa' run; if specified, parameters loaded from --load-params-file take priority over --s-value, --l-value, --pi-value, --sig2-zeroA, --sig2-zeroL values; note that --go-all-genes-label is used to populate initial sig2_gene values for all genes except base")
        parser.add_argument('--load-baseline-params-file', type=str, default=None, help="(optional) params of the baseline model to load and use for fold enrichment estimation; if --load-baseline-params-file is not provided, --load-params-file will be used;")

        parser.add_argument("--calc-detailed-trajectory", default=False, action="store_true", help=argparse.SUPPRESS) # "Calcualte detailed trajectory (one data point for each chromosome and --adam-epoch). Make cause excessively large .json files with --go-annot.")
        parser.add_argument("--calc-trajectory", default=False, action="store_true", help=argparse.SUPPRESS)
        parser.add_argument("--power-curve", default=False, action="store_true", help=argparse.SUPPRESS) # "Calculate power curve")
        parser.add_argument("--qqplot", default=False, action="store_true", help=argparse.SUPPRESS) #"Calculate QQ plots")
        parser.add_argument("--make-snps-file", default=False, action="store_true", help="a flag allowing to generate file with per-SNP estimates; will generate <out>.snps.csv output file")
        parser.add_argument('--downsample-factor', default=50, type=int, help=argparse.SUPPRESS)         # Applies to --power-curve and --qq-plots, --downsample-factor N' imply that only 1 out of N available z-score values will be used in model calculations.
        parser.add_argument('--se-samples', default=100, type=int, help="number of samples in standard error computation (using log-likelihood hessian / observed Fisher's information; applies only to --go-file-test")
        parser.add_argument("--gsa-base", default=False, action="store_true", help="set default parameters to values recommended for the baseline GSA model")
        parser.add_argument("--gsa-full", default=False, action="store_true", help="set default parameters to values recommended for the full GSA model")

    if analysis_type in ['fit1', 'test1', 'fit2', 'test2']:
        parser.add_argument("--annot-file", type=str, default=None, help="(optional) path to binary annotations in LD score regression format, i.e. <path>/baseline.@.annot.gz for fitting enrichment model model. This must include the first column with all ones ('base' annotation category covering the entire genome).")
        parser.add_argument("--go-file", type=str, default=None, help="(optional) path to GO antology file for fitting enrichment model model. The format is described in the documentation. 'base' category that covers entire genome will be added automatically.")
        parser.add_argument("--go-extend-bp", type=int, default=GO_EXTEND_BP_DEFAULT, help="extends each gene by this many base pairs, defining a symmetric window up and downstream (default: %(default)s)")
        parser.add_argument('--nckoef', default=None, type=float, help="(optional) coefficient translating total number of causal variants to the number of causal variants explaining 90%% of trait's heritability; by default this is estimated from the data (up until MiXeR v1.3 this was set to 0.319)")

        if analysis_type in ['fit1', 'test1']:
            parser.add_argument('--cubature-rel-error', type=float, default=1e-5, help="relative error for cubature stop criteria (default: %(default)s); applies to 'convolve' cost calculator")
            parser.add_argument('--cubature-max-evals', type=float, default=1000, help="max evaluations for cubature stop criteria (default: %(default)s); applies to 'convolve' cost calculator. "
                "Bivariate cubature require in the order of 10^4 evaluations and thus is much slower than sampling, therefore it is not exposed via mixer.py command-line interface. ")

        if analysis_type == 'fit1':
            parser.add_argument('--analysis', type=str, default='fit1', help=argparse.SUPPRESS, choices=['fit1', 'fit2', 'test1', 'test2'])
            parser.add_argument('--fit-sequence', type=str, default=['diffevo-fast', 'neldermead'], nargs='+', help=argparse.SUPPRESS, choices=['diffevo', 'diffevo-fast', 'neldermead', 'neldermead-fast', 'infinitesimal'])
            parser.add_argument('--load-params-file', type=str, default=None, help="params of the fitted model, for example from 'plsa' model (optional argument)")
        elif analysis_type == 'fit2':
            parser.add_argument('--analysis', type=str, default='fit2', help=argparse.SUPPRESS, choices=['fit1', 'fit2', 'test1', 'test2'])
            parser.add_argument('--fit-sequence', type=str, default=['diffevo-fast', 'neldermead-fast', 'brute1', 'brent1'], nargs='+', help=argparse.SUPPRESS, choices=['diffevo-fast', 'neldermead-fast', 'brute1', 'brute1-fast', 'brent1', 'brent1-fast', 'infinitesimal'])
            parser.add_argument('--trait1-params-file', type=str, default=None, help="univariate params for the first trait (required argument); applies to 'fit2' analysis only")
            parser.add_argument('--trait2-params-file', type=str, default=None, help="univariate params for the second trait (required argument); applies to 'fit2' analysis only")
        elif analysis_type == 'test1':
            parser.add_argument('--analysis', type=str, default='test1', help=argparse.SUPPRESS, choices=['fit1', 'fit2', 'test1', 'test2'])
            parser.add_argument('--fit-sequence', type=str, default=['load', 'inflation'], nargs='+', help=argparse.SUPPRESS, choices=['load', 'inflation'])
        elif analysis_type == 'test2':
            parser.add_argument('--analysis', type=str, default='test2', help=argparse.SUPPRESS, choices=['fit1', 'fit2', 'test1', 'test2'])
            parser.add_argument('--fit-sequence', type=str, default=['load', 'inflation'], nargs='+', help=argparse.SUPPRESS, choices=['load', 'inflation'])
        else:
            raise ValueError('internal error: invalid combination of do_fit and num_traits')

        if analysis_type in ['test1', 'test2']:
            parser.add_argument('--load-params-file', type=str, default=None, help="params of the fitted model (required argument); applies to 'test1' and 'test2' analyses")

        # all arguments below marked with argparse.SUPRESS option are internal and not recommended for a general use.
        # Valid options for --fit-sequence (remember to combine them in a sequence that makes sence):
        #       'load' reads previosly fitted parameters from a file (--load-params-file);
        #       'diffevo' performs a single iteration of differential evolution, which is the right way to setup an initial approximation;
        #       'neldermead' applies Nelder-Mead downhill simplex search;
        #       'brute1' applies only to bivariate optimization; it performs brute-force for one-dimentional optimization, searching optimal pi12 value constrained on genetic correlation (rg) and intercept (rho0);
        #       'brent1' is similar to brute1, but it uses brent method (inverse parabolic interpolation);
        #       'inflation' fits sig2zero (univariate) and rho_zero (bivariate), using fast cost function; this is quite special optimization step, typically useful for adjusting inflation parameters to another reference;
        #       'infinitesimal' fits a model with pi1=1 (univariate) or pi12=1 (bivariate) constrains, using fast cost function; this is quite special optimization step, typically used internally for AIC/BIC computation;
        # Note that bivariate fit is always constrained on univariate parameters, except for 'inflation' fit which adjust rho_zero and sig2_zeroA.
        # The '...-fast' optimizations use fast cost function.
        # Note that univariate optimization uses 'convolve' cost calculator, bivariate optimization uses 'sampling' cost calculator.
        # Typical univariate sequence: 'diffevo-fast neldermead'
        # Typical bivariate sequence: 'diffevo neldermead brute1 brent1'

        parser.add_argument('--diffevo-fast-repeats', type=int, default=20, help=argparse.SUPPRESS)      # repeat --diffevo-fast step this many times and choose the best run
        parser.add_argument('--downsample-factor', default=50, type=int, help=argparse.SUPPRESS)         # Applies to --power-curve and --qq-plots, --downsample-factor N' imply that only 1 out of N available z-score values will be used in model calculations.

        do_fit = (analysis_type in ['fit1', 'fit2'])
        parser.add_argument("--make-snps-file", default=False, action="store_true", help="Generate file with per-SNP estimates")

        parser.add_argument('--power-curve', default=(not do_fit), action="store_true", help=argparse.SUPPRESS) # generate power curves
        parser.add_argument('--qq-plots', default=(not do_fit), action="store_true", help=argparse.SUPPRESS)    # generate qq plot curves

        # Replace with The Sandwich Estimator ? http://www.stat.umn.edu/geyer/5601/notes/sand.pdf and re-implement CI for bivariate analysis?
        if analysis_type in ['fit1', 'test1']:
            parser.add_argument('--ci-alpha', type=float, default=None, help=argparse.SUPPRESS)              # significance level for the confidence interval estimation
            parser.add_argument('--ci-samples', type=int, default=10000, help=argparse.SUPPRESS)             # number of samples in uncertainty estimation
            parser.add_argument('--ci-power-samples', type=int, default=100, help=argparse.SUPPRESS)         # number of samples in power curves uncertainty estimation

    if analysis_type == 'snps':
        parser.add_argument('--r2', type=float, default=None, help="LD r2 threshold for prunning SNPs (default: %(default)s)")
        parser.add_argument('--maf', type=float, default=None, help="minor allele frequence (MAF) threshold (default: %(default)s)")
        parser.add_argument('--subset', type=int, default=2000000, help="number of SNPs to randomly select (default: %(default)s)")
        parser.add_argument('--seed', type=np.uint32, default=None, help="random seed (default: %(default)s)")
        parser.set_defaults(func=func)

    parser.set_defaults(func=func)

def parser_ld_add_arguments(func, parser):
    parser.add_argument("--bfile", type=str, default=None, help="path to plink bfile (required argument)")
    parser.add_argument('--r2min', type=float, default=0.05, help="r2 values above this threshold will be stored in sparse LD format (default: %(default)s)")
    parser.add_argument('--ldscore-r2min', type=float, default=0.001, help="r2 values above this threshold (and below --r2min) will be stored as LD scores that contribute to the cost function via an infinitesimal model (default: %(default)s)")
    parser.add_argument('--ld-window-kb', type=float, default=0, help="limit window similar to --ld-window-kb in 'plink r2'; 0 will disable this constraint (default: %(default)s); either ld-window-kb or ld-window argument must be provided")
    parser.add_argument('--ld-window', type=int, default=0, help="limit window similar to --ld-window in 'plink r2'; 0 will disable this constraint (default: %(default)s); either ld-window-kb or ld-window argument must be provided")
    parser.set_defaults(func=func)

def initialize_libbgmg_logging(args):
    libbgmg = LibBgmg(args.lib)
    log_fname = args.log if args.log else args.out + '.log'
    if os.path.exists(log_fname): os.system(f'rm {log_fname}')
    libbgmg.init_log(log_fname)
    log_header(args, sys.argv[1], libbgmg)
    libbgmg.dispose()
    return libbgmg

def generate_args_parser(__version__=None):
    parser = argparse.ArgumentParser(description=f"MiXeR v{__version__}: Univariate and Bivariate Causal Mixture for GWAS.")
    parser.add_argument('--version', action='version', version=f'MiXeR v{__version__}')

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--argsfile', type=open, action=LoadFromFile, default=None, help=argparse.SUPPRESS) #"file with additional command-line arguments, e.g. those that are across all of your mixer.py runs (--lib, --bim-file and --ld-file)")
    parent_parser.add_argument("--out", type=str, default="mixer", help="prefix for the output files, such as <out>.json (default: %(default)s);")
    parent_parser.add_argument("--lib", type=str, default="libbgmg.so", help="path to libbgmg.so plugin (default: %(default)s); can be also specified via BGMG_SHARED_LIBRARY env variable.")
    parent_parser.add_argument("--log", type=str, default=None, help="file to output log (default: <out>.log); NB! if --log points to an existing file the new lines will be appended to it at the end of the file.")

    subparsers = parser.add_subparsers()

    parser_ld_add_arguments(func=execute_ld_parser, parser=subparsers.add_parser("ld", parents=[parent_parser], help='prepare files with linkage disequilibrium information'))
    parser_add_arguments(func=execute_snps_parser, parser=subparsers.add_parser("snps", parents=[parent_parser], help='generate random sets of SNPs'), analysis_type='snps')
    parser_add_arguments(func=execute_fit1_or_test1_parser, parser=subparsers.add_parser("fit1", parents=[parent_parser], help='fit univariate MiXeR model'), analysis_type="fit1")
    parser_add_arguments(func=execute_fit1_or_test1_parser, parser=subparsers.add_parser("test1", parents=[parent_parser], help='test univariate MiXeR model'), analysis_type="test1")
    parser_add_arguments(func=execute_fit2_or_test2_parser, parser=subparsers.add_parser("fit2", parents=[parent_parser], help='fit bivariate MiXeR model'), analysis_type="fit2")
    parser_add_arguments(func=execute_fit2_or_test2_parser, parser=subparsers.add_parser("test2", parents=[parent_parser], help='test bivariate MiXeR model'), analysis_type="test2")
    parser_add_arguments(func=execute_plsa_parser, parser=subparsers.add_parser("plsa", parents=[parent_parser], help='fit PLSA model'), analysis_type='plsa')
    parser_add_arguments(func=execute_perf_parser, parser=subparsers.add_parser("perf", parents=[parent_parser], help='run performance evaluation of the MiXeR'), analysis_type='perf')
    parser_add_arguments(func=execute_fixed_effects_parser, parser=subparsers.add_parser("fixed_effects", parents=[parent_parser], help='find fixed effects (beta) at the level of tagging variants'), analysis_type='fixed_effects')
    parser_add_arguments(func=execute_split_sumstats_parser, parser=subparsers.add_parser("split_sumstats", parents=[parent_parser], help='split summary statistics file per-chromosome'), analysis_type='split_sumstats')

    return parser

def log_header(args, subparser_name, lib):
    defaults = vars(generate_args_parser().parse_args([subparser_name]))
    opts = vars(args)
    non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
    header = MASTHEAD
    header += "Call: \n"
    header += './mixer.py {} \\\n'.format(subparser_name)
    options = ['\t--'+x.replace('_','-')+' '+(' '.join([str(y) for y in opts[x]]) if isinstance(opts[x], list) else str(opts[x])).replace('\t', '\\t')+' \\' for x in non_defaults]
    header += '\n'.join(options).replace('True','').replace('False','')
    header = header[0:-1]+'\n'
    lib.log_message(header)

def load_params_file(fname_or_params, libbgmg, args):
    if isinstance(fname_or_params, AnnotUnivariateParams) or isinstance(fname_or_params, BivariateParams):
        return fname_or_params  # used for unit tests
    if isinstance(fname_or_params, str):
        data = json.loads(open(fname_or_params).read())
        return _dict_to_params(data['params'], libbgmg, args)
    else:
        return _dict_to_params(fname_or_params, libbgmg, args)

def apply_adam(parametrization, libbgmg_vec, args, vec, results, results_label, json_dump=False):
    adam_step_vec = []
    for epoch, step in zip(args.adam_epoch, args.adam_step): adam_step_vec.extend([step] * epoch)
    num_epoch=len(adam_step_vec)
    beta_1 = args.adam_beta1; beta_2 = args.adam_beta2; epsilon = args.adam_eps
    m = 0; v = 0; vec_step = None

    trajectory_detail = {'cost':[], 'vec':[], 'chr':[], 'step_size':[], 'm':[], 'v':[], 'm_hat':[], 'v_hat':[], 'gradients':[], 'params':[] }
    trajectory_fit = {'cost': [], 'params':[], 'step_size':[], 'annot_h2':[], 'annot_snps':[], 'annot_enrich':[], 'gene_h2':[], 'gene_snps':[],  'gene_enrich':[], 'total_h2':[], 'total_snps':[] }
    if args.calc_detailed_trajectory: results[results_label + '_detail'] = trajectory_detail
    if args.calc_trajectory: results[results_label + '_fit'] = trajectory_fit
    for t, step_size in zip(range(1, (num_epoch + 1)), adam_step_vec):
        fit_cost = 0; annot_h2 = 0; annot_snps = 0; gene_h2 = 0; gene_snps = 0
        for chri_index, chri in enumerate(np.random.permutation(len(libbgmg_vec))): 

            if vec_step is None: vec_step = [step_size for x in vec]

            libbgmg = libbgmg_vec[chri]
            parametrization._constraint._libbgmg = libbgmg

            pi_param = parametrization.vec_to_params(vec).pi

            libbgmg.set_option('kmax', 1 if (pi_param == 1) else args.kmax[0])
            libbgmg.set_option('seed', np.random.randint(np.iinfo(np.int32).max))
            cost, gradients = parametrization.calc_cost_with_gradients(vec, vec_step, verbose=False)

            m = beta_1 * m + (1 - beta_1) * gradients
            v = beta_2 * v + (1 - beta_2) * np.power(gradients, 2)
            m_hat = m / (1 - np.power(beta_1, t)) + (1 - beta_1) * gradients / (1 - np.power(beta_1, t))
            v_hat = v / (1 - np.power(beta_2, t))
            vec_next = list((np.array(vec).reshape([-1, 1]) + step_size * np.divide(m_hat, (np.sqrt(v_hat) + epsilon))).flatten())
            vec_step = list(np.array(vec_next) - np.array(vec))
            vec = vec_next

            params = parametrization.vec_to_params(vec)
            libbgmg.log_message('{}-{} {}'.format(t, chri_index + 1, params.as_string()))

            fit_cost += cost
            annot_h2 = annot_h2 + params.find_annot_h2(libbgmg.annomat).flatten()
            annot_snps = annot_snps + params.find_annot_snps(libbgmg.annomat).flatten()
            gene_h2 = gene_h2 + params.find_annot_h2(libbgmg.genemat).flatten()
            gene_snps = gene_snps + params.find_annot_snps(libbgmg.genemat).flatten()

            if args.calc_detailed_trajectory:
                trajectory_detail['cost'].append(cost)
                trajectory_detail['vec'].append(vec)
                trajectory_detail['chr'].append(str(chri))
                trajectory_detail['step_size'].append(step_size)
                trajectory_detail['m'].append(m)
                trajectory_detail['v'].append(v)
                trajectory_detail['m_hat'].append(m_hat)
                trajectory_detail['v_hat'].append(v_hat)
                trajectory_detail['gradients'].append(gradients)
                trajectory_detail['params'].append(_params_to_dict(params))

        if args.calc_trajectory:
            total_h2 = annot_h2.flat[0]; total_snps = annot_snps.flat[0]
            trajectory_fit['cost'].append(fit_cost)
            trajectory_fit['params'].append(_params_to_dict(params))
            trajectory_fit['step_size'].append(step_size)
            trajectory_fit['total_h2'].append(total_h2)
            trajectory_fit['total_snps'].append(total_snps)
            trajectory_fit['annot_h2'].append(annot_h2)
            trajectory_fit['annot_snps'].append(annot_snps)
            trajectory_fit['annot_enrich'].append(np.divide(np.divide(annot_h2, total_h2), np.divide(annot_snps, total_snps)))
            trajectory_fit['gene_h2'].append(gene_h2)
            trajectory_fit['gene_snps'].append(gene_snps)
            trajectory_fit['gene_enrich'].append(np.divide(np.divide(gene_h2, total_h2), np.divide(gene_snps, total_snps)))

        if json_dump:
            with open(args.out + '.json', 'w') as outfile:
                json.dump(results, outfile, cls=NumpyEncoder)

    return vec

def apply_plsa_univariate_fit_sequence(args, params, results, libbgmg_vec, trait_index):
    libbgmg = libbgmg_vec[0]

    if params is None:
        # historically, the initial value for sig2_annot was set to 1e-07, while for sig2_gene it was set to 1.0
        # note that for an arbitrary value "kappa>0" models (kappa*sig2_annot, sig2_gene/kappa) are equivalent.
        # now change initial values sot hat both sig2gene and sig2_annot are set to sqrt(1e-07) = 0.000316
        # the point is to have reasonably similar regularization coefficients for both sig2_annot and sig2_gene
        sig2_init_const = 0.000316

        params = AnnotUnivariateParams(pi=args.pi_value, s=args.s_value, l=args.l_value,
                                    sig2_zeroA=args.sig2_zeroA, sig2_zeroL=args.sig2_zeroL,
                                    sig2_beta=1, 
                                    sig2_annot=[sig2_init_const for i in libbgmg.annonames], p_annot=args.annot_p,
                                    sig2_gene=[sig2_init_const for i in libbgmg.genenames], p_gene=args.gene_p,
                                    libbgmg=libbgmg)

    # calculate h2 calibration constant, to set initial values of the sig2_gene or sig2_annot to a reasonable level
    if args.h2_init_calibration is not None:
        h2_calibration = np.sum([params.copy(libbgmg=libbgmg).find_h2() for libbgmg in libbgmg_vec])
        params._sig2_beta = params._sig2_beta * (args.h2_init_calibration / h2_calibration)
        libbgmg.log_message(f'h2_calibration={h2_calibration:.2e}')

    sig2_annot_const = [None for i in libbgmg.annonames] if ('annot' in args.fit) else params._sig2_annot
    sig2_gene_const = [None for i in libbgmg.genenames] if ('gene' in args.fit) else params._sig2_gene
    constraint = AnnotUnivariateParams(pi=None if ('pi' in args.fit) else params._pi,
                                       s=None if ('s' in args.fit) else params._s,
                                       l=None if ('l' in args.fit) else params._l,
                                       sig2_zeroA=None if ('sig2-zeroA' in args.fit) else params._sig2_zeroA,
                                       sig2_zeroL=None if ('sig2-zeroL' in args.fit) else params._sig2_zeroL,
                                       sig2_beta=None if (('annot' not in args.fit)) and (('gene' not in args.fit)) else params._sig2_beta,
                                       sig2_annot=sig2_annot_const,
                                       sig2_gene=sig2_gene_const, 
                                       p_annot=args.annot_p,
                                       p_gene=args.gene_p,
                                       libbgmg = libbgmg)
    parametrization = AnnotUnivariateParametrization(trait_index, constraint)
    if args.constrain_base_gene_category:
        parametrization._sig2_gene_base = params._sig2_gene[0]
    if args.nullify_go_all_genes_sig2_gene:
        parametrization._nullify_indices = [idx for idx, gene in enumerate(libbgmg.genenames) if (gene==args.go_all_genes_label)]

    libbgmg.log_message('apply --fit ' + ', '.join(args.fit))
    vec = parametrization.params_to_vec(params)
    if args.adam_separate_by_chromosome:
        for libbgmg in libbgmg_vec:
            chr_label = libbgmg._context_id
            libbgmg.log_message(f'--adam-separate-by-chromosome, chr {chr_label}')
            assert (args.z1max is None)
            if args.enable_faster_gradient_calculation: parametrization.enable_faster_gradient_computation(params, libbgmg, chr_label, args.go_all_genes_label)
            vec = apply_adam(parametrization, [libbgmg], args, vec, results, 'optimize_trajectory_inf_chr{chr_label}', json_dump=False)
            if args.enable_faster_gradient_calculation: parametrization.disable_faster_gradient_computation()
        params = parametrization.vec_to_params(vec)
        params._sig2_gene = [max(s, params._sig2_gene[0]) for s in params._sig2_gene]  # give second chance to parameters that otherwise are likely to be flushed to zero
        vec = parametrization.params_to_vec(params)
        for libbgmg in libbgmg_vec:
            chr_label = libbgmg._context_id
            libbgmg.log_message(f'--adam-separate-by-chromosome, chr {chr_label}')
            assert (args.z1max is None)
            if args.enable_faster_gradient_calculation: parametrization.enable_faster_gradient_computation(params, libbgmg, chr_label, args.go_all_genes_label)
            vec = apply_adam(parametrization, [libbgmg], args, vec, results, 'optimize_trajectory_inf2_chr{chr_label}', json_dump=False)
            if args.enable_faster_gradient_calculation: parametrization.disable_faster_gradient_computation()
    else:
        vec = apply_adam(parametrization, libbgmg_vec, args, vec, results, 'optimize_trajectory_inf', json_dump=False)
        params = parametrization.vec_to_params(vec)
        params._sig2_gene = [max(s, params._sig2_gene[0]) for s in params._sig2_gene]  # give second chance to parameters that otherwise are likely to be flushed to zero
        vec = parametrization.params_to_vec(params)
        vec = apply_adam(parametrization, libbgmg_vec, args, vec, results, 'optimize_trajectory_inf2', json_dump=False)
    params = parametrization.vec_to_params(vec)

    #for libbgmg in libbgmg_vec:
    #    parametrization._constraint._libbgmg = libbgmg
    #    cost, gradients = parametrization.calc_cost_with_gradients(vec, None, verbose=False)
    #    print(f'chr_label={libbgmg._context_id} cost={cost} gradients={gradients.flatten()}')

    return params

def apply_univariate_fit_sequence(args, libbgmg, fit_sequence, init_params=None, trait=1):
    params=init_params if (init_params is not None) else AnnotUnivariateParams(pi=1.0, libbgmg=libbgmg, s=0, l=0, sig2_zeroL=0)
    optimize_result_sequence=[]
    for fit_type in fit_sequence:
        libbgmg.log_message("fit_type=={}...".format(fit_type))

        if fit_type == 'load':
            libbgmg.log_message("Loading {}...".format(args.load_params_file))
            params = load_params_file(args.load_params_file, libbgmg, args)
            optimize_result = {}
            libbgmg.log_message("fit_type==load: Done, {}".format(params))

        elif (fit_type == 'diffevo') or fit_type == ('diffevo-fast'):
            libbgmg.set_option('cost_calculator', _cost_calculator_convolve if (fit_type == 'diffevo') else _cost_calculator_gaussian)
            parametrization = UnivariateParametrization_natural_axis(params=params, lib=libbgmg, trait=trait)

            h2_calibration = params.copy(pi=1.0, sig2_beta=1.0, sig2_zeroA=1.0).find_h2()
            h2_min = 1e-4; h2_max = 1
            sig2_beta_min = h2_min / h2_calibration
            sig2_beta_max = h2_max / h2_calibration
            totalhet = float(2.0 * np.dot(libbgmg._mafvec, 1.0 - libbgmg._mafvec))
            libbgmg.log_message(f'h2_calibration={h2_calibration:.2e}; totalhet={totalhet:.2e}; h2_min={h2_min:.2e}; h2_max={h2_max:.2e}; sig2_beta_min={sig2_beta_min:.2e}; sig2_beta_max={sig2_beta_max:.2e}')

            bounds_left = parametrization.params_to_vec(params.copy(pi=5e-5, sig2_beta=sig2_beta_min, sig2_zeroA=0.9))
            bounds_right = parametrization.params_to_vec(params.copy(pi=5e-1, sig2_beta=sig2_beta_max, sig2_zeroA=2.5))
            bounds4opt = [(l, r) for l, r in zip(bounds_left, bounds_right)]
            repeats = (1 if (fit_type == 'diffevo') else args.diffevo_fast_repeats)
            for repeat in range(repeats):
                optimize_result_tmp = scipy.optimize.differential_evolution(lambda x: parametrization.calc_cost(x), bounds4opt,
                    tol=0.01, mutation=(0.5, 1), recombination=0.7, seed=(args.seed + repeat), atol=0, updating='immediate', polish=False, workers=1)  #, **global_opt_options)
                params_tmp = parametrization.vec_to_params(optimize_result_tmp.x)
                libbgmg.log_message("--diffevo-fast-repeat={}: {} (cost={})".format(repeat, params_tmp, optimize_result_tmp.fun))
                if (repeat == 0) or (optimize_result_tmp.fun < optimize_result.fun):
                    optimize_result = optimize_result_tmp
                    selected_repeat = repeat
            libbgmg.log_message("--diffevo-fast-repeat={} has the lowest cost function and will be used".format(selected_repeat))
            params = parametrization.vec_to_params(optimize_result.x).copy()

        elif (fit_type == 'neldermead') or (fit_type == 'neldermead-fast'):
            if params == None: raise(RuntimeError('params == None, unable to proceed apply "neldermead" fit'))
            libbgmg.set_option('cost_calculator', _cost_calculator_convolve if (fit_type == 'neldermead') else _cost_calculator_gaussian)
            parametrization = UnivariateParametrization_natural_axis(params=params, lib=libbgmg, trait=trait)
            optimize_result = scipy.optimize.minimize(lambda x: parametrization.calc_cost(x), parametrization.params_to_vec(params),
                method='Nelder-Mead', options={'maxiter':1200, 'fatol':1e-7, 'xatol':1e-4, 'adaptive':True})
            params = parametrization.vec_to_params(optimize_result.x).copy()

        elif fit_type == 'inflation':
            if params == None: raise(RuntimeError('params == None, unable to proceed apply "inflation" fit'))
            libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
            parametrization = UnivariateParametrization_constPI_constSIG2BETA(params=params, lib=libbgmg, trait=trait)
            optimize_result = scipy.optimize.minimize(lambda x: parametrization.calc_cost(x), parametrization.params_to_vec(params), method='Nelder-Mead')
            params = parametrization.vec_to_params(optimize_result.x).copy()

        elif fit_type == 'infinitesimal':
            if params == None: raise(RuntimeError('params == None, unable to proceed apply "infinitesimal" fit'))
            libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
            inft_params = params.copy(pi=1.0, sig2_beta=params._sig2_beta)
            parametrization = UnivariateParametrization_constPI(params=inft_params, lib=libbgmg, trait=trait)
            optimize_result = scipy.optimize.minimize(lambda x: parametrization.calc_cost(x), parametrization.params_to_vec(inft_params), method='Nelder-Mead')
            params = parametrization.vec_to_params(optimize_result.x).copy()

        else:
            libbgmg.log_message("Unable to apply {} in univariate fit".format(fit_type))

        libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
        if optimize_result:
            enhance_optimize_result(optimize_result, cost_n=np.sum(libbgmg.weights), cost_fast=params.cost(libbgmg, trait))
        optimize_result['params']=_params_to_dict(params)   # params after optimization
        optimize_result_sequence.append((fit_type, optimize_result))
        libbgmg.log_message("fit_type=={} done ({}, {})".format(fit_type, params, optimize_result))

    if params == None: raise(RuntimeError('Empty --fit-sequence'))
    return params, optimize_result_sequence

def apply_bivariate_fit_sequence(args, libbgmg, fit_sequence):
    if args.analysis == 'fit2':
        libbgmg.log_message("Loading univariate constrains for bivariate analysis...")
        params1 = load_params_file(args.trait1_params_file, libbgmg, args)
        params2 = load_params_file(args.trait2_params_file, libbgmg, args)
        libbgmg.log_message("trait1: {}".format(params1))
        libbgmg.log_message("trait2: {}".format(params2))

    params = None; optimize_result_sequence = []
    for fit_type in fit_sequence:
        libbgmg.log_message("fit_type=={}...".format(fit_type))

        if fit_type == 'load':
            libbgmg.log_message("Loading {}...".format(args.load_params_file))
            params = load_params_file(args.load_params_file, libbgmg, args)
            params1 = params._params1
            params2 = params._params2
            libbgmg.log_message("fit_type==load... trait1 params: {}".format(params1))            
            libbgmg.log_message("fit_type==load... trait2 params: {}".format(params2))    
            libbgmg.log_message("fit_type=={} done ({})".format(fit_type, params))
            continue

        elif fit_type == 'inflation':
            if params == None: raise(RuntimeError('params == None, unable to proceed with "{}" fit'.format(fit_type)))
            params1, ors1 = apply_univariate_fit_sequence(args, libbgmg, ['inflation'], init_params=params1, trait=1)
            optimize_result_sequence.extend(ors1)
            params2, ors2 = apply_univariate_fit_sequence(args, libbgmg, ['inflation'], init_params=params2, trait=2)
            optimize_result_sequence.extend(ors2)
            libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
            parametrization = BivariateParametrization_constUNIVARIATE_constRHOBETA_constPI(
                const_params1=params1, const_params2=params2, const_pi12=params._pi12, const_rho_beta=params._rho_beta, lib=libbgmg)
            optimize_result = scipy.optimize.minimize(lambda x: parametrization.calc_cost(x), parametrization.params_to_vec(params), method='Nelder-Mead')
            params = parametrization.vec_to_params(optimize_result.x)

        elif (fit_type == 'diffevo') or fit_type == ('diffevo-fast'):
            libbgmg.set_option('cost_calculator', _cost_calculator_sampling if (fit_type == 'diffevo') else _cost_calculator_gaussian)
            parametrization = BivariateParametrization_constUNIVARIATE_natural_axis(lib=libbgmg, const_params1=params1, const_params2=params2)
            max_pi12 = min(params1.pi, params2.pi)
            bounds_left = parametrization.params_to_vec(BivariateParams(params1=params1, params2=params2, pi12=0.05 * max_pi12, rho_beta=-0.95, rho_zero=-0.95))
            bounds_right = parametrization.params_to_vec(BivariateParams(params1=params1, params2=params2, pi12=0.95 * max_pi12, rho_beta=0.95, rho_zero=0.95))
            bounds4opt = [(l, r) for l, r in zip(bounds_left, bounds_right)]
            repeats = (1 if (fit_type == 'diffevo') else args.diffevo_fast_repeats)
            for repeat in range(repeats):
                optimize_result_tmp = scipy.optimize.differential_evolution(lambda x: parametrization.calc_cost(x), bounds4opt,
                    tol=0.01, mutation=(0.5, 1), recombination=0.7, seed=(args.seed + repeat), atol=0, updating='immediate', polish=False, workers=1)  #, **global_opt_options)
                params_tmp = parametrization.vec_to_params(optimize_result_tmp.x)
                libbgmg.log_message("--diffevo-fast-repeat={}: {} (cost={})".format(repeat, params_tmp, optimize_result_tmp.fun))
                if (repeat == 0) or (optimize_result_tmp.fun < optimize_result.fun):
                    optimize_result = optimize_result_tmp
                    selected_repeat = repeat
            libbgmg.log_message("--diffevo-fast-repeat={} has the lowest cost function and will be used".format(selected_repeat))
            params = parametrization.vec_to_params(optimize_result.x)

        elif (fit_type == 'neldermead') or (fit_type == 'neldermead-fast'):
            if params == None: raise(RuntimeError('params == None, unable to proceed with "{}" fit'.format(fit_type)))
            libbgmg.set_option('cost_calculator', _cost_calculator_sampling if (fit_type == 'neldermead') else _cost_calculator_gaussian)
            parametrization = BivariateParametrization_constUNIVARIATE_natural_axis(lib=libbgmg, const_params1=params1, const_params2=params2)
            optimize_result = scipy.optimize.minimize(lambda x: parametrization.calc_cost(x), parametrization.params_to_vec(params),
                method='Nelder-Mead', options={'maxiter':1200, 'fatol':1e-7, 'xatol':1e-4, 'adaptive':True})
            params = parametrization.vec_to_params(optimize_result.x)

        elif (fit_type == 'brute1') or (fit_type == 'brute1-fast'):
            # brute1 optimization intentionally forgets previously fitted value of params._pi12, to avoid being stuck in a local minimum.
            if params == None: raise(RuntimeError('params == None, unable to proceed with "{}" fit'.format(fit_type)))
            libbgmg.set_option('cost_calculator', _cost_calculator_sampling if (fit_type == 'brute1') else _cost_calculator_gaussian)
            parametrization = BivariateParametrization_constUNIVARIATE_constRG_constRHOZERO(
                lib=libbgmg, const_params1=params1, const_params2=params2, const_rg=params._rg(), const_rho_zero=params._rho_zero)
            min_pi12 = parametrization._min_pi12; max_pi12 = parametrization._max_pi12
            params._pi12 = 0.99 * min_pi12 + 0.01 * max_pi12; bounds_left = parametrization.params_to_vec(params)
            params._pi12 = 0.01 * min_pi12 + 0.99 * max_pi12; bounds_right = parametrization.params_to_vec(params)
            bounds4opt = [(l, r) for l, r in zip(bounds_left, bounds_right)]
            x0, fval, grid, Jout = scipy.optimize.brute(lambda x: parametrization.calc_cost(x), bounds4opt, Ns=20, full_output=True)
            optimize_result = scipy.optimize.OptimizeResult(fun=fval, nit=1, nfev=len(grid),
                            success=True, message="Optimization terminated successfully.", x=x0)
            optimize_result['grid'] = grid
            optimize_result['Jout'] = Jout
            optimize_result['grid_params'] = [_params_to_dict(parametrization.vec_to_params([x])) for x in grid]
            params = parametrization.vec_to_params(optimize_result.x)

        elif (fit_type == 'brent1') or (fit_type == 'brent1-fast'):
            if params == None: raise(RuntimeError('params == None, unable to proceed with "{}" fit'.format(fit_type)))
            libbgmg.set_option('cost_calculator', _cost_calculator_sampling if (fit_type == 'brent1') else _cost_calculator_gaussian)
            parametrization = BivariateParametrization_constUNIVARIATE_constRG_constRHOZERO(
                lib=libbgmg, const_params1=params1, const_params2=params2, const_rg=params._rg(), const_rho_zero=params._rho_zero)
            min_pi12 = parametrization._min_pi12; max_pi12 = parametrization._max_pi12
            bracket_middle = parametrization.params_to_vec(params)
            params._pi12 = min_pi12; bounds_left = parametrization.params_to_vec(params)
            params._pi12 = max_pi12; bounds_right = parametrization.params_to_vec(params)
            bracket4opt = (bounds_left[0], bracket_middle[0], bounds_right[0])
            bracketfunc = (parametrization.calc_cost(bounds_left[0]), parametrization.calc_cost(bracket_middle[0]),parametrization.calc_cost(bounds_right[0]))
            libbgmg.log_message("bracket4opt = {}, {}".format(bracket4opt,bracketfunc))
            try:
                xmin, fval, iter, funcalls = scipy.optimize.brent(lambda x: parametrization.calc_cost(x), brack=bracket4opt, full_output=True)
                success=True; message = "Optimization terminated successfully."
            except ValueError as ve:
                xmin = bracket_middle[0]; fval = parametrization.calc_cost(bracket_middle[0])
                iter = 1; funcalls = 4 # three byscipy.optimize.brent, plus one in the line above
                success=False; message=str(ve)
            optimize_result = scipy.optimize.OptimizeResult(fun=fval, nit=iter, nfev=funcalls, x=[xmin], success=success, message=message)
            params = parametrization.vec_to_params(optimize_result.x)

        elif (fit_type == 'infinitesimal'):
            libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
            if (params1._pi != 1.0): libbgmg.log_message("change params1._pi from {} to 1.0".format(params1._pi))
            if (params2._pi != 1.0): libbgmg.log_message("change params2._pi from {} to 1.0".format(params2._pi))
            params1._pi = 1.0
            params2._pi = 1.0
            params = BivariateParams(params1=params1, params2=params2, pi12=1, rho_beta=0, rho_zero=0)
            parametrization = BivariateParametrization_constUNIVARIATE_infinitesimal(lib=libbgmg, const_params1=params1, const_params2=params2)
            optimize_result = scipy.optimize.minimize(lambda x: parametrization.calc_cost(x), parametrization.params_to_vec(params),
                method='Nelder-Mead', options={'maxiter':1200, 'fatol':1e-7, 'xatol':1e-4, 'adaptive':True})
            params = parametrization.vec_to_params(optimize_result.x)

        else:
            libbgmg.log_message("Unable to apply {} in bivariate fit".format(fit_type))

        libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
        # cost_df=9 --- nine free parameters (incl. univariate)
        enhance_optimize_result(optimize_result, cost_df=9, cost_n=np.sum(libbgmg.weights), cost_fast=params.cost(libbgmg))
        optimize_result['params']=_params_to_dict(params)   # params after optimization
        optimize_result_sequence.append((fit_type, optimize_result))
        libbgmg.log_message("fit_type=={} done ({}, {})".format(fit_type, params, optimize_result))

    if params == None: raise(RuntimeError('Empty --fit-sequence'))
    return (params, params1, params2, optimize_result_sequence)

# helper function to debug non-json searizable types...
def print_types(results, libbgmg):
    if isinstance(results, dict):
        for k, v in results.items():
            libbgmg.log_message('{}: {}'.format(k, type(v)))
            print_types(v, libbgmg)

def execute_ld_parser(args):
    initialize_libbgmg_logging(args)
    libbgmg = LibBgmg(args.lib)
    libbgmg.calc_ld_matrix(args.bfile, args.out, args.r2min, args.ldscore_r2min, args.ld_window, args.ld_window_kb)
    libbgmg.log_message('Done')

def find_annomat(annot_file, chr_labels, lib_file_suffix, libbgmg):
    if (annot_file != None):
        libbgmg.log_message('Loading annotations from {}...'.format(annot_file.replace('@', lib_file_suffix)))
        df = pd.concat([pd.read_csv(annot_file.replace('@', str(chr_label)), sep='\t') for chr_label in chr_labels])
        bim = df[[col for col in ['CHR', 'BP', 'SNP', 'CM'] if (col in df.columns)]]
        for col in ['CHR', 'BP', 'SNP', 'CM']:
            if col in df.columns:
                del df[col]

        annomat = csr_matrix(df.values.astype(np.float32))
        annonames = df.columns.values
        del df
        libbgmg.log_message('Done, {} annotations on {} SNPs available'.format(len(annonames), annomat.shape[0]))
    else:
        annomat = csr_matrix(np.ones(shape=(libbgmg.num_snp, 1), dtype=np.float32))
        annonames = ['base']
        bim = pd.DataFrame({'BP':np.zeros(shape=(libbgmg.num_snp,))})
    return (annomat, annonames, bim)

def load_go_file(go_file, go_extend_bp, exclude_genes=[]):
    go = pd.read_csv(go_file, sep='\t')
    if 'base' in go['GO'].values:
        raise ValueError(f'{go_file} file has "base" in its "GO" column; this is not allowed as "base" category is added automatically')

    go = go[~go['GENE'].isin(exclude_genes)].reset_index(drop=True)
    go_clean = go['CHR FROM TO GO'.split()].copy()
    if not np.all(go_clean['GO'].notnull()): raise(ValueError('encounter null values in "GO" column'))
    go_map = {val:(index+1) for index, val  in enumerate(go_clean['GO'].unique())}; go_map['base'] = 0
    go_clean['GO_index'] = go_clean['GO'].map(go_map)

    go_clean.loc[:, 'FROM'] = go_clean['FROM'].values + go_clean['CHR'].values * CHR_OFFSET - go_extend_bp
    go_clean.loc[:, 'TO'] = go_clean['TO'].values + go_clean['CHR'].values * CHR_OFFSET + go_extend_bp + 1   # +1 to include right end of the intervals
    go_clean = pd.concat([go_clean, pd.DataFrame({'CHR':[0], 'FROM':[0], 'TO':[100*CHR_OFFSET], 'GO':['base'], 'GO_index':[0]})])

    go_map_inv = {v: k for k, v in go_map.items()}
    genenames = [go_map_inv[x] for x in range(len(go_map_inv))]
    assert(len(genenames) == len(go_map))

    return go_clean, genenames

def find_genemat(go_file, go_extend_bp, bim_file, chr_labels, lib_file_suffix, libbgmg):
    if (go_file != None):
        libbgmg.log_message('Loading {} file(s)...'.format(bim_file.replace('@', lib_file_suffix)))
        bim = pd.concat([pd.read_csv(bim_file.replace('@', str(chr_label)), sep='\t', header=None, names='CHR SNP GP BP A1 A2'.split()) for chr_label in chr_labels])
        bim.reset_index(inplace=True, drop=True)
        snp_coordinate = bim['BP'].values + bim['CHR'].values * CHR_OFFSET

        libbgmg.log_message('Loading {} with --go-extend-bp={}...'.format(go_file, go_extend_bp))
        go_clean, genenames = load_go_file(go_file, go_extend_bp)

        libbgmg.log_message('Annotating SNPs to GO gene-sets...')
        coo = libbgmg.annotate_intervals(snp_coordinate,go_clean['FROM'].values,go_clean['TO'].values,go_clean['GO_index'].values)
        genemat = coo_matrix((np.ones(len(coo[0, :]), dtype=np.float32), (coo[0, :], coo[1, :])), shape=(len(bim), len(genenames))).tocsr()
        genemat[:, 0] = 1-(np.sum(genemat[:, 1:], 1) > 0).astype(int).flatten()
        genemat.eliminate_zeros()

        libbgmg.log_message('Done, annotation matrix has {} non-zero entries across {} SNP and {} GO terms'.format(np.int64(genemat.sum()), genemat.shape[0], genemat.shape[1]))
    else:
        genemat = csr_matrix(np.ones(shape=(libbgmg.num_snp, 1), dtype=np.float32))
        genenames = ['base']
        bim = None
    return (genemat, genenames, bim)

def initialize_mixer_plugin(args, context_id, trait1_file, trait2_file, randprune_n, loadlib_file, savelib_file, chr_labels, lib_file_suffix):
    exclude_ranges = ' '.join([f'{range.chr}:{int(range.from_bp)}-{int(range.to_bp)}' for range in make_ranges(args.exclude_ranges)]) if ('exclude_ranges' in args) else ""

    libbgmg = LibBgmgAnnot(args.lib, context_id=context_id)
    if loadlib_file:
        libbgmg.log_message(f'loading {loadlib_file}...')
        if args.ld_file: libbgmg.log_message(f'WARNING: --ld-file argument will be ignored')
        libbgmg.load(loadlib_file)

        if set(libbgmg.chrnumvec).symmetric_difference(set(chr_labels)):
            raise(ValueError(f'--chr-labels argument appear to be different than chromosomes {list(sorted(set(libbgmg.chrnumvec)))} loaded from --loadlib-file'))

        if 'allow_ambiguous_snps' in args:
            libbgmg.set_option('allow_ambiguous_snps', args.allow_ambiguous_snps)
        if trait1_file: libbgmg.read_trait_file(1, trait1_file, args.exclude if ('exclude' in args) else "", args.extract if ('extract' in args) else "", exclude_ranges)
        if trait2_file: libbgmg.read_trait_file(2, trait2_file, args.exclude if ('exclude' in args) else "", args.extract if ('extract' in args) else "", exclude_ranges)
        libbgmg.set_option('cost_calculator', _cost_calculator_sampling)

        for opt, val in convert_args_to_libbgmg_options(args):
            libbgmg.set_option(opt, val)

        set_weights(args, libbgmg, randprune_n)
        libbgmg.set_option('diag', 0)
    else:
        if 'use_complete_tag_indices' in args:
            libbgmg.set_option('use_complete_tag_indices', args.use_complete_tag_indices)
        if 'allow_ambiguous_snps' in args:
            libbgmg.set_option('allow_ambiguous_snps', args.allow_ambiguous_snps)
        libbgmg.init(args.bim_file, "", chr_labels,
                    trait1_file,
                    trait2_file,
                    args.exclude if ('exclude' in args) else "",
                    args.extract if ('extract' in args) else "",
                    exclude_ranges)
        libbgmg.set_option('cost_calculator', _cost_calculator_sampling)

        for opt, val in convert_args_to_libbgmg_options(args):
            libbgmg.set_option(opt, val)

        for chr_label in chr_labels: 
            libbgmg.set_ld_r2_coo_from_file(chr_label, args.ld_file.replace('@', str(chr_label)))
            libbgmg.set_ld_r2_csr(chr_label)

        if trait1_file or trait2_file: set_weights(args, libbgmg, randprune_n)
        libbgmg.set_option('diag', 0)
        if savelib_file:
            if not args.use_complete_tag_indices: libbgmg.log_message('WARNING: unless you know exactly what you are doing, it is recommended to use --use-complete-tag-indices together with --savelib-file.')
            libbgmg.log_message(f'saving {savelib_file}...')
            libbgmg.save(savelib_file)

    annomat, annonames, bim = find_annomat(args.annot_file if ('annot_file' in args) else None, chr_labels, lib_file_suffix, libbgmg)
    libbgmg.annomat = annomat
    libbgmg.annonames = annonames
    if bim is not None: libbgmg.bim = bim

    genemat, genenames, bim = find_genemat(args.go_file if ('go_file' in args) else None, args.go_extend_bp if ('go_extend_bp' in args) else None, args.bim_file, chr_labels, lib_file_suffix, libbgmg)
    libbgmg.genemat = genemat
    libbgmg.genenames = genenames
    if bim is not None: libbgmg.bim = bim

    return libbgmg

def execute_perf_parser(args):
    initialize_libbgmg_logging(args)
    fix_and_validate_args(args)
    libbgmg = initialize_mixer_plugin(args, 0, args.trait1_file, args.trait2_file, args.randprune_n, args.loadlib_file, args.savelib_file, args.chr2use, '@')
    params1 = UnivariateParams_obsolete(pi=3e-3, sig2_beta=1e-5, sig2_zeroA=1.05)
    params12 = BivariateParams(params1=params1, params2=params1, pi12=2e-3, rho_beta=0.5, rho_zero=0.3)

    perf_data = []
    for threads in args.threads:
        for kmax in args.kmax:
            for costcalc in _cost_calculator:
                libbgmg.set_option('threads', threads)
                libbgmg.set_option('kmax', kmax)
                libbgmg.set_option('cost_calculator', _cost_calculator[costcalc])

                if (kmax!=np.max(args.kmax)) and (costcalc in ['gaussian', 'convolve']):
                    continue  # gaussian and convolve are calculated only for one value of kmax, e.g. for the largest one
                if (costcalc == 'convolve') and (args.trait2_file):
                    continue  # skip convolve for bivariate analysis

                start = time.time()
                cost = (params12.cost(libbgmg) if args.trait2_file else params1.cost(libbgmg, trait=1))
                end = time.time()
                perf_data.append((threads, kmax, costcalc, end-start, cost))
                libbgmg.log_message('threads={}, kmax={}, costcalt={} took {} seconds'.format(threads, kmax, costcalc, end-start))

    pd.DataFrame(perf_data, columns=['threads', 'kmax', 'costcalc', 'time_sec', 'cost']).to_csv(args.out + '.csv', sep='\t', index=False)
    libbgmg.log_message('Done')

def execute_snps_parser(args):
    initialize_libbgmg_logging(args)
    fix_and_validate_args(args)

    libbgmg = initialize_mixer_plugin(args, 0, "", "", None, None, None, args.chr2use, '@')
    mafvec = np.minimum(libbgmg.mafvec, 1-libbgmg.mafvec)

    libbgmg.log_message('Load {}...'.format(args.bim_file))
    ref=pd.concat([pd.read_csv(args.bim_file.replace('@', str(chr_label)), sep='\t', header=None, names='CHR SNP GP BP A1 A2'.split()) for chr_label in args.chr2use])
    libbgmg.log_message('{} SNPs in total'.format(len(ref)))

    # step0 - generate random values for clumping (this is a way to implement random pruning)
    buf = np.random.rand(libbgmg.num_tag,)
    
    # step1 - filter SNPs below MAF threshold
    if args.maf:
        buf[mafvec<args.maf] = np.nan
        libbgmg.log_message('{} SNPs pass --maf {} threshold'.format(np.sum(np.isfinite(buf)), args.maf))
    
    # step2 - select a random subset of SNPs
    if args.subset > 0:
        indices = list(np.where(np.isfinite(buf))[0])
        sample = random.sample(indices, min(args.subset, len(indices)))
        mask = np.ones(buf.shape); mask[sample] = 0
        buf[mask == 1] = np.nan
        libbgmg.log_message('{} SNPs randomly selected (--subset)'.format(np.sum(np.isfinite(buf))))

    # step3 - further prune SNPs at certain r2 threshold
    if args.r2:
        buf = libbgmg.perform_ld_clump(args.r2, buf.flatten()) 
        libbgmg.log_message('{} SNPs pass random prunning at --r2 {} threshold'.format(np.sum(np.isfinite(buf)), args.r2))

    # step4 - save results
    ref.SNP[np.isfinite(buf)].to_csv(args.out, index=False, header=False)
    libbgmg.log_message('Result saved to {}'.format(args.out))

    # step 5 - save .info file
    #pd.DataFrame({'chrnumvec': libbgmg.chrnumvec, 'mafvec': libbgmg.mafvec, 'tldvec': libbgmg.ld_sum_r2}).to_csv(args.out + '.info', sep='\t', index=False)
    
    libbgmg.log_message('Done')

def init_results_struct(args, libbgmg_vec=None):
    results = {}
    results['options'] = vars(args).copy()
    if libbgmg_vec:
        results['options']['totalhet'] = np.sum([float(2.0 * np.dot(libbgmg._mafvec, 1.0 - libbgmg._mafvec)) for libbgmg in libbgmg_vec])
        results['options']['num_snp'] = np.sum([float(libbgmg.num_snp) for libbgmg in libbgmg_vec])
        results['options']['num_tag'] = np.sum([float(libbgmg.num_tag) for libbgmg in libbgmg_vec])
        results['options']['sum_weights'] = np.sum([float(np.sum(libbgmg.weights)) for libbgmg in libbgmg_vec])
        results['options']['trait1_nval'] = float(np.nanmedian(np.concatenate([libbgmg.get_nvec(trait=1) for libbgmg in libbgmg_vec])))
        results['options']['time_started'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        results['ci'] = {}
    return results

def execute_fit1_or_test1_parser(args):
    initialize_libbgmg_logging(args)
    fix_and_validate_args(args)

    libbgmg = initialize_mixer_plugin(args, 0, args.trait1_file, "", args.randprune_n, args.loadlib_file, args.savelib_file, args.chr2use, '@')

    if not args.trait1_file:
        libbgmg.log_message('--trait1-file not provided, exit')
        return

    df_weights = pd.DataFrame({'chrnumvec':libbgmg.chrnumvec, 'weights': tagvec_as_snpvec(libbgmg.weights, libbgmg.defvec)})
    results = init_results_struct(args, [libbgmg])
    results['analysis'] = 'univariate'

    if (args.load_params_file is not None) and (args.fit_sequence[0] != 'load'):
        args.fit_sequence = ['load'] + args.fit_sequence

    libbgmg.log_message('--fit-sequence: {}...'.format(args.fit_sequence))

    params, optimize_result = apply_univariate_fit_sequence(args, libbgmg, args.fit_sequence)
    results['params'] = _params_to_dict(params)
    results['optimize'] = optimize_result

    libbgmg.log_message('Calculate AIC/BIC w.r.t. infinitesimal model (fast cost function)...')
    params_inft, optimize_result_inft = apply_univariate_fit_sequence(args, libbgmg, ['infinitesimal'], init_params=params)
    results['inft_params'] = _params_to_dict(params_inft)
    results['inft_optimize'] = optimize_result_inft

    libbgmg.set_option('kmax', 1 if (params._pi == 1) else args.kmax[0])
    libbgmg.set_option('kmax_pdf', 1 if (params._pi == 1) else args.kmax_pdf[0])

    NCKoef = params.find_nckoef() if (args.nckoef is None) else args.nckoef
    results['options']['nckoef'] = NCKoef
    funcs, _ = _calculate_univariate_uncertainty_funcs(None, NCKoef)
    for func_name, func in funcs: results['ci'][func_name] = {'point_estimate': func(params)}

    ci_sample = None
    if args.ci_alpha and np.isfinite(args.ci_alpha):
        libbgmg.log_message("Uncertainty estimation...")
        libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
        parametrization = UnivariateParametrization_natural_axis(params, libbgmg, trait=1)
        results['ci'], ci_sample = _calculate_univariate_uncertainty(params, parametrization, args.ci_alpha, args.ci_samples, NCKoef)
        for k, v in results['ci'].items():
            libbgmg.log_message('{}: point_estimate={:.3g}, mean={:.3g}, median={:.3g}, std={:.3g}, ci=[{:.3g}, {:.3g}]'.format(k, v['point_estimate'], v['mean'], v['median'], v['std'], v['lower'], v['upper']))
        libbgmg.log_message("Uncertainty estimation done.")

    if args.power_curve:
        trait_index = 1
        results['power'] = calc_power_curve(libbgmg, params, trait_index, args.downsample_factor)
        if ci_sample is not None:
            power_ci = []
            for ci_index, ci_params in enumerate(ci_sample[:(args.ci_power_samples-1)]):
                libbgmg.log_message("Power curves uncertainty, {} of {}".format(ci_index, args.ci_power_samples))
                power_ci.append(calc_power_curve(libbgmg, ci_params, trait_index, args.downsample_factor))
            results['power_ci'] = power_ci

    if args.qq_plots:
        trait_index = 1
        mask = np.ones((libbgmg.num_tag, ), dtype=bool)

        mafvec = libbgmg._mafvec[libbgmg.defvec]
        hetvec = np.float32(2.0) * mafvec * (1-mafvec)
        tldvec = libbgmg._tldvec[libbgmg.defvec]

        tag_defvec = (libbgmg.weights > 0) & np.isfinite(libbgmg.zvec1) & np.isfinite(libbgmg.nvec1)

        results['qqplot'] = calc_qq_plot(libbgmg, params, trait_index, args.downsample_factor, mask,
            title='het \\in [{:.3g},{:.3g}); L \\in [{:.3g},{:.3g})'.format(np.min(hetvec),np.max(hetvec),np.min(tldvec),np.max(tldvec)))

        het_bins = np.concatenate(([-np.inf], np.quantile(hetvec[tag_defvec], [1/3, 2/3]), [np.inf]))
        tld_bins = np.concatenate(([-np.inf], np.quantile(tldvec[tag_defvec], [1/3, 2/3]), [np.inf]))
        het_bins_display = np.concatenate(([np.min(hetvec[tag_defvec])], np.quantile(hetvec[tag_defvec], [1/3, 2/3]), [np.max(hetvec[tag_defvec])]))
        tld_bins_display = np.concatenate(([np.min(tldvec[tag_defvec])], np.quantile(tldvec[tag_defvec], [1/3, 2/3]), [np.max(tldvec[tag_defvec])]))
        results['qqplot_bins'] = []
        for i in range(0, 3):
            for j in range(0, 3):
                mask = ((hetvec>=het_bins[i]) & (hetvec<het_bins[i+1]) & (tldvec >= tld_bins[j]) &  (tldvec < tld_bins[j+1]))
                results['qqplot_bins'].append(calc_qq_plot(libbgmg, params, trait_index, args.downsample_factor, mask,
                    title='het \\in [{:.3g},{:.3g}); L \\in [{:.3g},{:.3g})'.format(het_bins_display[i], het_bins_display[i+1], tld_bins_display[j], tld_bins_display[j+1])))

    if args.make_snps_file:
        trait_index = 1
        libbgmg.log_message('Finding per-SNP information...')
        df=find_per_snp_information_univariate([libbgmg], params, trait_index, args.make_snps_file)
        fname = args.out + '.snps.csv'
        df.to_csv(fname, sep='\t', index=False)
        libbgmg.log_message('--make-snps-file output saved to {}, covering {} SNPs'.format(fname, len(df)))

    results['options']['time_finished'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(args.out + '.json', 'w') as outfile:
        json.dump(results, outfile, cls=NumpyEncoder)
    if args.save_weights: df_weights.to_csv(args.out + '.weights', index=False, na_rep=0, sep='\t')

    libbgmg.set_option('diag', 0)
    libbgmg.log_message('Done')

def execute_fit2_or_test2_parser(args):
    initialize_libbgmg_logging(args)
    fix_and_validate_args(args)

    libbgmg = initialize_mixer_plugin(args, 0, args.trait1_file, args.trait2_file, args.randprune_n, args.loadlib_file, args.savelib_file, args.chr2use, '@')

    df_weights = pd.DataFrame({'chrnumvec':libbgmg.chrnumvec, 'weights': tagvec_as_snpvec(libbgmg.weights, libbgmg.defvec)})
    results = init_results_struct(args, [libbgmg])
    results['analysis'] = 'bivariate'
    results['options']['trait2_nval'] = float(np.nanmedian(libbgmg.get_nvec(trait=2)))

    libbgmg.log_message('--fit-sequence: {}...'.format(args.fit_sequence))
    params, _, _, optimize_result = apply_bivariate_fit_sequence(args, libbgmg, args.fit_sequence)
    results['params'] = _params_to_dict(params)
    results['optimize'] = optimize_result

    libbgmg.set_option('kmax', 1 if (params._pi12 == 1) else args.kmax[0])
    libbgmg.set_option('kmax_pdf', 1 if (params._pi12 == 1) else args.kmax_pdf[0])

    NCKoef = params.find_nckoef() if (args.nckoef is None) else args.nckoef
    results['options']['nckoef'] = NCKoef
    funcs, _ = _calculate_bivariate_uncertainty_funcs(None, NCKoef)
    for func_name, func in funcs: results['ci'][func_name] = {'point_estimate': func(params)}

    if args.qq_plots:
        zgrid, pdf = calc_bivariate_pdf(libbgmg, params, args.downsample_factor)
        results['pdf_zgrid'] = zgrid
        results['pdf'] = pdf
        results['qqplot'] = calc_bivariate_qq(libbgmg, zgrid, pdf)

    if args.make_snps_file:
        libbgmg.log_message('Finding per-SNP information...')
        df=find_per_snp_information_bivariate(libbgmg, params)
        fname = args.out + '.snps.csv'
        df.to_csv(fname, sep='\t', index=False)
        libbgmg.log_message('--make-snps-file output saved to {}, covering {} SNPs'.format(fname, len(df)))

    results['options']['time_finished'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(args.out + '.json', 'w') as outfile:
        json.dump(results, outfile, cls=NumpyEncoder)
    if args.save_weights: df_weights.to_csv(args.out + '.weights', index=False, na_rep=0, sep='\t')

    libbgmg.set_option('diag', 0)
    libbgmg.log_message('Done')

def tagvec_as_snpvec(tagvec, defvec):
        snpvec = np.empty(defvec.shape, dtype=float)
        snpvec[:] = np.nan
        snpvec[defvec] = tagvec
        return snpvec

def find_per_snp_information_bivariate(libbgmg, params):
    # Save a file with the following information
    original_weights = libbgmg.weights
    posterior_weights = np.ones(original_weights.shape, dtype=np.float32)
    for trait_index in [1, 2]:
        posterior_weights[~np.isfinite(libbgmg.get_zvec(trait_index)) | ~np.isfinite(libbgmg.get_nvec(trait_index))] = 0
    libbgmg.weights = posterior_weights
    c00, c10, c01, c20, c11, c02 = params.delta_posterior(libbgmg)

    # find tag_pdf (full model) as well as 5 reduced models
    # - associations with shared component
    # - associatinos with first component
    # - associatinos with second component
    # - associations with both unique components, but not through the shared component
    # - null associations
    num_snp = libbgmg.num_snp
    pi_mat_orig = params.find_pi_mat(num_snp)

    libbgmg.set_option('aux_option', _auxoption_tagpdf)
    libbgmg.set_option('cost_calculator', _cost_calculator_sampling)

    tag_pdf_xxx = libbgmg.calc_unified_bivariate_aux(pi_mat_orig, params.find_sig2_mat(), params.find_rho_vec(num_snp),
        sig2_zeroA=[params._params1._sig2_zeroA, params._params2._sig2_zeroA],
        sig2_zeroC=[1, 1], sig2_zeroL=[0, 0], rho_zeroA=params._rho_zero, rho_zeroL=0).flat[:libbgmg.num_tag]

    pi_mat = pi_mat_orig.copy(); pi_mat[:, :] = 1e-10
    tag_pdf_ooo = libbgmg.calc_unified_bivariate_aux(pi_mat, params.find_sig2_mat(), params.find_rho_vec(num_snp),
        sig2_zeroA=[params._params1._sig2_zeroA, params._params2._sig2_zeroA],
        sig2_zeroC=[1, 1], sig2_zeroL=[0, 0], rho_zeroA=params._rho_zero, rho_zeroL=0).flat[:libbgmg.num_tag]

    pi_mat = pi_mat_orig.copy(); pi_mat[1, :] = 1e-10; pi_mat[2, :] = 1e-10
    tag_pdf_xoo = libbgmg.calc_unified_bivariate_aux(pi_mat, params.find_sig2_mat(), params.find_rho_vec(num_snp),
        sig2_zeroA=[params._params1._sig2_zeroA, params._params2._sig2_zeroA],
        sig2_zeroC=[1, 1], sig2_zeroL=[0, 0], rho_zeroA=params._rho_zero, rho_zeroL=0).flat[:libbgmg.num_tag]

    pi_mat = pi_mat_orig.copy(); pi_mat[0, :] = 1e-10; pi_mat[2, :] = 1e-10
    tag_pdf_oxo = libbgmg.calc_unified_bivariate_aux(pi_mat, params.find_sig2_mat(), params.find_rho_vec(num_snp),
        sig2_zeroA=[params._params1._sig2_zeroA, params._params2._sig2_zeroA],
        sig2_zeroC=[1, 1], sig2_zeroL=[0, 0], rho_zeroA=params._rho_zero, rho_zeroL=0).flat[:libbgmg.num_tag]

    pi_mat = pi_mat_orig.copy(); pi_mat[0, :] = 1e-10; pi_mat[1, :] = 1e-10
    tag_pdf_oox = libbgmg.calc_unified_bivariate_aux(pi_mat, params.find_sig2_mat(), params.find_rho_vec(num_snp),
        sig2_zeroA=[params._params1._sig2_zeroA, params._params2._sig2_zeroA],
        sig2_zeroC=[1, 1], sig2_zeroL=[0, 0], rho_zeroA=params._rho_zero, rho_zeroL=0).flat[:libbgmg.num_tag]

    pi_mat = pi_mat_orig.copy(); pi_mat[0, :] = 1e-10
    tag_pdf_oxx = libbgmg.calc_unified_bivariate_aux(pi_mat, params.find_sig2_mat(), params.find_rho_vec(num_snp),
        sig2_zeroA=[params._params1._sig2_zeroA, params._params2._sig2_zeroA],
        sig2_zeroC=[1, 1], sig2_zeroL=[0, 0], rho_zeroA=params._rho_zero, rho_zeroL=0).flat[:libbgmg.num_tag]

    pi_mat = pi_mat_orig.copy(); pi_mat[1, :] = 1e-10
    tag_pdf_xox = libbgmg.calc_unified_bivariate_aux(pi_mat, params.find_sig2_mat(), params.find_rho_vec(num_snp),
        sig2_zeroA=[params._params1._sig2_zeroA, params._params2._sig2_zeroA],
        sig2_zeroC=[1, 1], sig2_zeroL=[0, 0], rho_zeroA=params._rho_zero, rho_zeroL=0).flat[:libbgmg.num_tag]

    pi_mat = pi_mat_orig.copy(); pi_mat[2, :] = 1e-10
    tag_pdf_xxo = libbgmg.calc_unified_bivariate_aux(pi_mat, params.find_sig2_mat(), params.find_rho_vec(num_snp),
        sig2_zeroA=[params._params1._sig2_zeroA, params._params2._sig2_zeroA],
        sig2_zeroC=[1, 1], sig2_zeroL=[0, 0], rho_zeroA=params._rho_zero, rho_zeroL=0).flat[:libbgmg.num_tag]

    libbgmg.set_option('aux_option', _auxoption_none)

    libbgmg.weights = original_weights
    defvec = libbgmg.defvec
    df_snps = pd.DataFrame({
        'defvec': defvec,
        'chrnumvec': libbgmg.chrnumvec,
        'sig2beta_T1': params._params1.find_sig2_mat().flatten(),
        'sig2beta_T2': params._params2.find_sig2_mat().flatten(),
        'mafvec': libbgmg._mafvec,
        'tldvec': libbgmg._tldvec,
        'weights': tagvec_as_snpvec(libbgmg.weights, defvec),
        'zvec1': tagvec_as_snpvec(libbgmg.get_zvec(trait=1), defvec),
        'nvec1': tagvec_as_snpvec(libbgmg.get_nvec(trait=1), defvec),
        'zvec2': tagvec_as_snpvec(libbgmg.get_zvec(trait=2), defvec),
        'nvec2': tagvec_as_snpvec(libbgmg.get_nvec(trait=2), defvec),
        'posteriordeltaC00': tagvec_as_snpvec(c00, defvec),
        'posteriordeltaC10': tagvec_as_snpvec(c10, defvec),
        'posteriordeltaC01': tagvec_as_snpvec(c01, defvec),
        'posteriordeltaC20': tagvec_as_snpvec(c20, defvec),
        'posteriordeltaC11': tagvec_as_snpvec(c11, defvec),
        'posteriordeltaC02': tagvec_as_snpvec(c02, defvec),
        'tag_pdf_xxx_full': tagvec_as_snpvec(tag_pdf_xxx, defvec),
        'tag_pdf_ooo_null': tagvec_as_snpvec(tag_pdf_ooo, defvec),
        'tag_pdf_xoo_uniqueT1': tagvec_as_snpvec(tag_pdf_xoo, defvec),
        'tag_pdf_oxo_uniqueT2': tagvec_as_snpvec(tag_pdf_oxo, defvec),
        'tag_pdf_oox_shared': tagvec_as_snpvec(tag_pdf_oox, defvec),
        'tag_pdf_xxo_twohits': tagvec_as_snpvec(tag_pdf_xxo, defvec),
        'tag_pdf_xox_notuniqueT2': tagvec_as_snpvec(tag_pdf_xox, defvec),
        'tag_pdf_oxx_notuniqueT1': tagvec_as_snpvec(tag_pdf_oxx, defvec),
    })
    return df_snps

def find_per_snp_information_univariate(libbgmg_vec, params, trait_index, extended=False):
    # Save a file with the following information
    df_snps = []
    for libbgmg in libbgmg_vec:
        params._libbgmg = libbgmg
        params._snp_defvec = None
        if extended:
            original_weights = libbgmg.weights
            posterior_weights = np.zeros(original_weights.shape, dtype=np.float32)
            posterior_weights[np.isfinite(libbgmg.get_zvec(trait_index)) & np.isfinite(libbgmg.get_nvec(trait_index))] = 1
            libbgmg.weights = posterior_weights
            c0, c1, c2 = params.delta_posterior(libbgmg, trait_index)
            libbgmg.set_option('cost_calculator', _cost_calculator_sampling)
            tag_pdf = params.tag_pdf(libbgmg, trait_index)
            tag_ez2 = params.tag_ez2(libbgmg, trait_index)
            libbgmg.weights = original_weights
            defvec = libbgmg.defvec
            df_snps.append(pd.DataFrame({
                'defvec': defvec,
                'chrnumvec': libbgmg.chrnumvec,
                'h2vec': params.find_h2_vec().flatten(),
                'sig2beta': params.find_sig2_mat().flatten(),
                'mafvec': libbgmg._mafvec,
                'tldvec': libbgmg._tldvec,
                'weights': tagvec_as_snpvec(libbgmg.weights, defvec),
                'zvec': tagvec_as_snpvec(libbgmg.get_zvec(trait_index), defvec),
                'nvec': tagvec_as_snpvec(libbgmg.get_nvec(trait_index), defvec),
                'posteriordeltaC0': tagvec_as_snpvec(c0, defvec),
                'posteriordeltaC1': tagvec_as_snpvec(c1, defvec),
                'posteriordeltaC2': tagvec_as_snpvec(c2, defvec),
                'tag_pdf': tagvec_as_snpvec(tag_pdf, defvec),
                'tag_ez2': tagvec_as_snpvec(tag_ez2, defvec),
            }))
        else:
            df_snps.append(pd.DataFrame({
                'chrnumvec': libbgmg.chrnumvec,
                'h2vec': params.find_h2_vec().flatten(),
                'sig2beta': params.find_sig2_mat().flatten(),
                'mafvec': libbgmg._mafvec,
            }))

    return pd.concat(df_snps)

def find_snp_defvec(chr_label, posvec, exclude_ranges, libbgmg):
    snp_defvec = None

    for range in exclude_ranges:
        if str(chr_label) != str(range.chr): continue
        snp_defvec2 = (posvec < range.from_bp) | (posvec > range.to_bp)
        if snp_defvec is None: snp_defvec = snp_defvec2
        else: snp_defvec = snp_defvec & snp_defvec2

    if snp_defvec is not None:
        snp_defvec = np.array(snp_defvec, dtype=np.float32).reshape((-1, 1))
        libbgmg.log_message("Exclude {} SNPs on chr{}".format(int(len(snp_defvec) - np.sum(snp_defvec)), chr_label))

    return snp_defvec

def compute_annot_enrichment(annot_file_test, chr2use, exclude_ranges, libbgmg_vec, params, base_params):
    libbgmg_vec[0].log_message("Compute annot enrichment for {}...".format(annot_file_test))
    annot_h2 = 0; annot_h2_base = 0; annot_snps = 0
    for chr_label, libbgmg in zip(chr2use, libbgmg_vec):
        annomat, annonames, bim = find_annomat(annot_file_test, [str(chr_label)], str(chr_label), libbgmg)

        snp_defvec = find_snp_defvec(chr_label, bim['BP'], exclude_ranges, libbgmg)
        for p in [params]:
            p._libbgmg = libbgmg
            p._snp_defvec = snp_defvec
        
        annot_h2 = annot_h2 + params.find_annot_h2(annomat).flatten()

        h2vec_base = base_params[base_params['chrnumvec'] == chr_label]['h2vec'].values.reshape((-1, 1))
        annot_h2_base = annot_h2_base + annomat.transpose().dot(h2vec_base).flatten()

        hetvec = np.float32(2.0) * params._mafvec * (1-params._mafvec)
        annot_snps = annot_snps + annomat.transpose().dot(hetvec)

    if 'base' not in annonames: raise(ValueError('"base" not found among annonames'))
    idx_all_genes = [i for i, g in enumerate(annonames) if g=='base'][0]
    total_h2 = annot_h2.flat[idx_all_genes]
    total_h2_base = annot_h2_base.flat[idx_all_genes]
    annot_enrich  = np.divide(np.divide(annot_h2, total_h2), np.divide(annot_h2_base, total_h2_base))
    df_enrich = pd.DataFrame({'h2':annot_h2.flat, 'h2_base':annot_h2_base.flat, 'enrich':annot_enrich.flat, 'snps':annot_snps.flat, 'annonames':annonames})
    return df_enrich

def compute_geneset_enrichment(go_file_test, args, exclude_ranges, libbgmg_vec, params_vec, base_params, calc_loglike_diff=False):
    libbgmg_vec[0].log_message("Compute geneset enrichment for {}...".format(go_file_test))

    # first element is point estimates; the rest are samples from posterior distribution to estimate error bars
    params = params_vec[0]

    gene_h2 = None; gene_h2_base = 0; gene_snps = 0
    loglike_diff = None
    for chr_label, libbgmg in zip(args.chr2use, libbgmg_vec):
        genemat, genenames, bim = find_genemat(go_file_test, args.go_extend_bp, args.bim_file, [str(chr_label)], str(chr_label), libbgmg)

        if gene_h2 is None: # initialize
            gene_h2 = np.zeros(shape=(len(params_vec), len(genenames)))

        snp_defvec = find_snp_defvec(chr_label, bim['BP'], exclude_ranges, libbgmg)
        for p in params_vec:
            p._libbgmg = libbgmg
            p._snp_defvec = snp_defvec

        for params_idx, p in enumerate(params_vec):
            gene_h2[params_idx, :] = gene_h2[params_idx, :] + p.find_annot_h2(genemat).flatten()

        h2vec_base = base_params[base_params['chrnumvec'] == chr_label]['h2vec'].values.reshape((-1, 1))
        gene_h2_base = gene_h2_base + genemat.transpose().dot(h2vec_base).flatten()
        hetvec = np.float32(2.0) * params._mafvec * (1-params._mafvec)
        gene_snps = gene_snps + genemat.transpose().dot(hetvec)

        # ====================== frankenstein cost start ======================
        if loglike_diff is None:
            loglike_diff = np.zeros((len(genenames), 1))

        if calc_loglike_diff:
            trait_index = 1
            libbgmg.set_option('kmax', 1 if (params._pi == 1) else args.kmax[0])
            pi_mat_full = params.find_pi_mat(libbgmg.num_snp)
            sig2_mat_full = params.find_sig2_mat()
            cost_full = libbgmg.calc_unified_univariate_cost(trait_index, pi_mat_full, sig2_mat_full, params._sig2_zeroA, sig2_zeroC=1, sig2_zeroL=params._sig2_zeroL)

            sig2_mat_base = base_params[base_params['chrnumvec'] == chr_label]['sig2beta'].values.reshape((-1, 1))

            for gene_idx, gene in enumerate(genenames):
                gene_snp_indices = genemat[:, gene_idx].nonzero()[0]
                if len(gene_snp_indices) == 0: continue

                sig2_mat_frank = sig2_mat_full.copy()
                sig2_mat_frank[gene_snp_indices] = sig2_mat_base[gene_snp_indices]

                pi_mat_frank = pi_mat_full
                #pi_mat_frank = pi_mat_full.copy()
                #pi_mat_frank[gene_snp_indices] = pi_mat_base[gene_snp_indices]

                cost_frank = libbgmg.calc_unified_univariate_cost(trait_index, pi_mat_frank, sig2_mat_frank, params._sig2_zeroA, sig2_zeroC=1, sig2_zeroL=params._sig2_zeroL)
                loglike_diff[gene_idx] += (cost_frank - cost_full)
        # ====================== frankenstein cost finish ======================

    go_all_genes_label = args.go_all_genes_label if (args.go_all_genes_label in genenames) else 'base'
    idx_all_genes = [i for i, g in enumerate(genenames) if g==go_all_genes_label][0]
    total_h2 = gene_h2[0, idx_all_genes]
    total_h2_base = gene_h2_base.flat[idx_all_genes]
    gene_enrich  = np.divide(np.divide(gene_h2[0, :], total_h2), np.divide(gene_h2_base, total_h2_base))
    gene_enrich_std  = np.divide(np.divide(np.std(gene_h2[1:, :], axis=0), total_h2), np.divide(gene_h2_base, total_h2_base))
    df_enrich = pd.DataFrame({'h2':gene_h2[0, :].flat, 'se_h2':np.std(gene_h2[1:, :], axis=0),  'h2_base':gene_h2_base.flat, 'enrich':gene_enrich.flat, 'se_enrich':gene_enrich_std, 'snps':gene_snps.flat, 'GO':genenames, 'loglike_diff':loglike_diff.flat})
    df_enrich = pd.merge(df_enrich, pd.read_csv(go_file_test, sep='\t')[['GO','GENE']].groupby('GO').count().reset_index(), on='GO', how='left').rename(columns={'GENE':'NGENES'})
    return df_enrich

def execute_plsa_parser(args):
    libbgmg = initialize_libbgmg_logging(args)
    fix_and_validate_args(args)
    if (args.load_params_file is None) and (args.adam_disable):
        raise ValueError('--load-params-file is required for --adam-disable')

    if args.go_all_genes_label == 'base': raise ValueError('"--go-all-genes-label base" is not allowed')

    if args.gsa_base:
        args.fit = ['sig2-zeroA', 'l', 'annot', 'gene']
        args.l_value = -0.25
        args.h2_init_calibration = 0.1
        args.se_samples = 0
        libbgmg.log_message(f'apply --gsa-base settings (--fit {" ".join(args.fit)} --l-value {args.l_value} --h2-init-calibration {args.h2_init_calibration})')

    if args.gsa_full:
        args.fit = ['gene']
        args.constrain_base_gene_category = True
        args.nullify_go_all_genes_sig2_gene = True
        args.calc_loglike_diff_go_test = 'fast'
        args.adam_separate_by_chromosome = True
        args.adam_beta1 = 0.1
        args.adam_beta2 = 0.8
        libbgmg.log_message(f'apply --gsa-full settings (--fit {" ".join(args.fit)}  --constrain-base-gene-category --nullify-go-all-genes-sig2-gene --calc-loglike-diff-go-test fast --adam-separate-by-chromosome --adam-beta1 {args.adam_beta1} --adam-beta2 {args.adam_beta2})')

    if args.power_curve:
        raise NotImplementedError('power curve not implemented for "mixer.py plsa". Use "mixer.py test1 --load-params-file <out>.json", where <out>.json is generated by "mixer.py plsa" call.')
    if args.save_ldsc_reference and args.trait1_file:
        # --trait1-file will be ignored due to --save-ldsc-reference flag
        args.trait1_file = ""

    if args.annot_file_test or args.go_file_test:
        if (not args.load_params_file) and (not args.load_baseline_params_file):
            raise ValueError('--annot-file-test/--go-file-test require either --load-params-file or --load-baseline-params-file to be provided')
        if (args.load_baseline_params_file is None) and (args.load_params_file is not None):
            args.load_baseline_params_file = args.load_params_file
        check_input_file(args, 'load_baseline_params_file')

    libbgmg_vec = []
    annot_snps = 0; gene_snps = 0
    for chr_label in args.chr2use:
        libbgmg = initialize_mixer_plugin(
            args,
            context_id=chr_label,
            trait1_file=args.trait1_file.replace('@', str(chr_label)),
            trait2_file=args.trait2_file.replace('@', str(chr_label)) if ('trait2_file' in args) else "",
            randprune_n=args.randprune_n, 
            loadlib_file=args.loadlib_file.replace('@', str(chr_label)) if args.loadlib_file else None, 
            savelib_file=args.savelib_file.replace('@', str(chr_label)) if args.savelib_file else None, 
            chr_labels=[chr_label], 
            lib_file_suffix=str(chr_label))
        annonames = libbgmg.annonames
        genenames = libbgmg.genenames
        libbgmg_vec.append(libbgmg)
        annot_snps = annot_snps + np.sum(libbgmg.annomat, 0).flatten()
        gene_snps = gene_snps + np.sum(libbgmg.genemat, 0).flatten()

        if args.save_ldsc_reference:
            savemat_dict = {}

            if args.annot_file and args.go_file:
                raise(ValueError('for --save-ldsc-reference, use either --annot-file or --go-file, but not both at the same time'))
            if len(genenames) > 2000:
                raise(ValueError('Too many genenames to save into LDSC format'))

            bim = pd.read_csv(args.bim_file.replace('@', str(chr_label)), sep='\t', header=None, names='CHR SNP GP BP A1 A2'.split())

            tag_array, snp_array, r_array = libbgmg.get_ld_r_chr(chr_label)
            r2=coo_matrix((np.power(r_array, 2), (tag_array, snp_array)))
            m_5_50_idx = (libbgmg._mafvec > 0.05) & (libbgmg._mafvec < 0.95)

            savemat_dict['ld_idx1'] = tag_array
            savemat_dict['ld_idx2'] = snp_array
            savemat_dict['ld_r'] = r_array
            savemat_dict['mafvec'] = libbgmg._mafvec
            savemat_dict['chrnumvec'] = bim['CHR'].values
            savemat_dict['posvec'] = bim['BP'].values
            savemat_dict['rsidvec'] = bim['SNP'].values
            savemat_dict['A1'] = bim['A1'].values
            savemat_dict['A2'] = bim['A2'].values

            if args.annot_file:
                fname = args.out.replace('@', str(chr_label)) + '.l2.ldscore.gz'
                libbgmg.log_message(f'saving {fname}...')
                annomat_r2 = np.dot(r2, libbgmg.annomat).todense()
                annomat_r2_df = pd.concat([libbgmg.bim[['CHR', 'SNP', 'BP']], 
                    pd.DataFrame(annomat_r2, columns=[x+'L2' for x in libbgmg.annonames])], axis=1)
                annomat_r2_df.to_csv(fname, sep='\t', index=False)

                savemat_dict['annonames'] = libbgmg.annonames
                savemat_dict['annomat'] = libbgmg.annomat.todense()
                savemat_dict['annomat_ld'] = annomat_r2

                open(args.out.replace('@', str(chr_label)) + '.l2.M', 'w').write(
                    '\t'.join([str(int(x)) for x in np.sum(libbgmg.annomat, 0).flat]))
                open(args.out.replace('@', str(chr_label)) + '.l2.M_5_50', 'w').write(
                    '\t'.join([str(int(x)) for x in np.sum(libbgmg.annomat[m_5_50_idx, :], 0).flat]))

            if args.go_file:
                fname = args.out.replace('@', str(chr_label)) + '.annot.gz'
                libbgmg.log_message(f'saving {fname}...')
                genemat_df = pd.concat([libbgmg.bim[['CHR', 'SNP', 'BP', 'GP']].rename(columns={'GP':'CM'}), 
                    pd.DataFrame(libbgmg.genemat.todense(), columns=libbgmg.genenames)], axis=1)
                del genemat_df['base']
                genemat_df.to_csv(fname, sep='\t', index=False)

                fname = args.out.replace('@', str(chr_label)) + '.l2.ldscore.gz'
                libbgmg.log_message(f'saving {fname}...')
                genemat_r2 = np.dot(r2, libbgmg.genemat).todense()
                genemat_r2_df = pd.concat([libbgmg.bim[['CHR', 'SNP', 'BP']], 
                    pd.DataFrame(genemat_r2, columns=[x+'L2' for x in libbgmg.genenames])], axis=1)
                del genemat_r2_df['baseL2']
                genemat_r2_df.to_csv(fname, sep='\t', index=False)

                open(args.out.replace('@', str(chr_label)) + '.l2.M', 'w').write(
                    '\t'.join([str(int(x)) for x in np.sum(libbgmg.genemat[:, 1:], 0).flat]))
                open(args.out.replace('@', str(chr_label)) + '.l2.M_5_50', 'w').write(
                    '\t'.join([str(int(x)) for x in np.sum(libbgmg.genemat[m_5_50_idx, 1:], 0).flat]))

            sio.savemat(args.out.replace('@', str(chr_label)) + '.mat', savemat_dict, format='5', do_compression=False, oned_as='column', appendmat=False)

    if not args.trait1_file:
        libbgmg.log_message('--trait1-file not provided, exit')
        return

    if 0 in annot_snps:
        without_snps = ', '.join([annoname for num_snps, annoname in zip(list(annot_snps.flat), annonames) if (num_snps == 0)])
        libbgmg.log_message('WARNING: the following entries from --annot-file contain 0 SNPs: {}'.format(without_snps))
    if 0 in gene_snps:
        without_snps = ', '.join([genename for num_snps, genename in zip(list(gene_snps.flat), genenames) if (num_snps == 0)])
        libbgmg.log_message('WARNING: the following entries from --go-file contain 0 SNPs: {}'.format(without_snps))

    trait_index = 1
    df_weights = pd.concat([pd.DataFrame({'chrnumvec':libbgmg.chrnumvec, 'weights': tagvec_as_snpvec(libbgmg.weights, libbgmg.defvec)}) for libbgmg in libbgmg_vec])
    results = init_results_struct(args, libbgmg_vec)
    results['analysis'] = 'plsa'
    results['options']['annonames'] = annonames
    results['options']['genenames'] = genenames
    if args.load_params_file is not None:
        libbgmg.log_message("Loading {}...".format(args.load_params_file))
        params = load_params_file(args.load_params_file, libbgmg, args)
        libbgmg.log_message(params.as_string())
    else:
        params = None

    if not args.adam_disable:
        params = apply_plsa_univariate_fit_sequence(args, params, results, libbgmg_vec, trait_index)
        results['options']['time_finished_fit'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        libbgmg.log_message(params.as_string())
    else:
        libbgmg.log_message("Skip optimiztion due to --adam-disable")

    results['params'] = _params_to_dict(params)

    if args.gsa_base or args.make_snps_file:
        libbgmg.log_message('Finding per-SNP information...')
        df=find_per_snp_information_univariate(libbgmg_vec, params, trait_index, args.make_snps_file)
        fname = args.out + '.snps.csv'
        df.to_csv(fname, sep='\t', index=False)
        libbgmg.log_message('--make-snps-file output saved to {}, covering {} SNPs'.format(fname, len(df)))

    if args.qqplot:
        libbgmg.log_message('Making QQ plot...')
        results['qqplot'] = calc_qq_plot_plsa(libbgmg_vec, params, trait_index, args.downsample_factor)
        fname = args.out + '.qq.json'
        with open(fname, 'w') as outfile:
            json.dump(results['qqplot'], outfile, cls=NumpyEncoder)
        libbgmg.log_message('--qqplot output saved to {}'.format(fname))

    if args.se_samples > 0: # use observed Fisher's information to calculate error bars
        sig2_gene_vec_samples = np.tile([x for x in params._sig2_gene], (args.se_samples, 1))
        sig2_gene_vec_samples_orig = sig2_gene_vec_samples.copy()

        for chr_label, libbgmg in zip(args.chr2use, libbgmg_vec):
            params._libbgmg = libbgmg
            params._snp_defvec = None

            B, gene_indices_on_chr = AnnotUnivariateParametrization.calc_faster_gradient_computation(
                libbgmg, params, chr_label, args.go_all_genes_label)
            libbgmg.log_message(f'Computing hessian for {len(gene_indices_on_chr)} genes on chr {chr_label}')

            sig2_gene_chr = params._sig2_gene[gene_indices_on_chr]
            sig2_vec = B.dot(sig2_gene_chr) + params.sig2_zeroA

            weights_hess = np.multiply(libbgmg.weights, np.divide(sig2_vec - 2 * np.power(libbgmg.zvec1, 2), 2 * np.power(sig2_vec, 3)))
            hess_chr = B.transpose().dot(scipy.sparse.diags(weights_hess, dtype=np.float32).dot(B))  # hessian w.r.t. sig2_gene

            neg_hess_chr = -hess_chr.todense()
            neg_hess_chr = (neg_hess_chr+neg_hess_chr.T)/2

            eig_vals, eig_vecs = np.linalg.eig(neg_hess_chr)
            eig_vals_sorted = np.array(list(sorted(eig_vals[eig_vals>0])))
            eig_vals_indices = np.where(np.cumsum(eig_vals_sorted) / np.sum(eig_vals_sorted) < 1e-4)[0]
            eig_value_thr = eig_vals_sorted[np.max(eig_vals_indices)] if len(eig_vals_indices > 0) else np.min(eig_vals_sorted)
            eig_vals_indices = np.where(np.cumsum(eig_vals_sorted) / np.sum(eig_vals_sorted) < 1e-8)[0]
            eig_value_min = eig_vals_sorted[np.max(eig_vals_indices)] if len(eig_vals_indices > 0) else np.min(eig_vals_sorted)

            for eig_idx in np.where(eig_vals < eig_value_thr)[0]:
                nonSPD_indices = np.where(np.abs(eig_vecs[:, eig_idx])>0.1)[0]
                for nonSPD_index in nonSPD_indices:
                    diag = neg_hess_chr[nonSPD_index, nonSPD_index]
                    neg_hess_chr[nonSPD_index, :] = 0
                    neg_hess_chr[:, nonSPD_index] = 0
                    neg_hess_chr[nonSPD_index, nonSPD_index] = diag
                    if neg_hess_chr[nonSPD_index, nonSPD_index] < eig_value_min:
                        neg_hess_chr[nonSPD_index, nonSPD_index] = eig_value_min

            neg_hess_chr_inv = np.linalg.inv(neg_hess_chr)
            sig2_gene_vec_samples[:, gene_indices_on_chr] += np.random.multivariate_normal(np.zeros(neg_hess_chr.shape[0]),
                                            neg_hess_chr_inv,
                                            size=args.se_samples, check_valid='warn', tol=1e-8)

        se_factor = 5  # as discussed in the final round of revisions
        sig2_gene_vec_samples = np.minimum(np.maximum(1e-15, sig2_gene_vec_samples), se_factor*sig2_gene_vec_samples_orig)
        sig2_gene_vec_samples = np.maximum(1e-15, sig2_gene_vec_samples)
        sig2_gene_samples = sig2_gene_vec_samples

    results['options']['time_finished'] = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(args.out + '.json', 'w') as outfile:
        json.dump(results, outfile, cls=NumpyEncoder)
    if args.save_weights: df_weights.to_csv(args.out + '.weights', index=False, na_rep=0, sep='\t')
    libbgmg.log_message('Done, results saved to {}'.format(args.out + '.json'))

    base_params = None
    if args.load_baseline_params_file is not None:
        fname = args.load_baseline_params_file.rstrip('.json') + '.snps.csv'
        libbgmg.log_message(f"Loading {fname} ...")
        base_params = pd.read_csv(fname, sep='\t')

    # calculate fold enrichments of heritability
    exclude_ranges = make_ranges(args.exclude_ranges)
    if args.annot_file_test:
        df_enrich = compute_annot_enrichment(args.annot_file_test, args.chr2use, exclude_ranges, libbgmg_vec, params, base_params)
        df_enrich.to_csv(args.out + '.annot_test_enrich.csv', sep='\t', index=False)
        libbgmg.log_message("Enrichment results for {} are saved to {}".format(args.annot_file_test, args.out + '.annot_test_enrich.csv'))

    if args.go_file_test:
        # prepare params_vec for SE computations
        params_vec = [params]
        for i in range(args.se_samples):
            p = params.copy()
            p._sig2_gene = sig2_gene_samples[i, :]
            params_vec.append(p)

        # compute gene-set enrichments, SE's, loglike_diff and loglike_df
        df_enrich_test = compute_geneset_enrichment(args.go_file_test, args, exclude_ranges, libbgmg_vec, params_vec, base_params, args.calc_loglike_diff_go_test == 'full')

        # for the --calc-loglike-diff-go-test 'fast' algorithm we compute enrichment for --go-file and then aggregate to gene-sets
        if args.calc_loglike_diff_go_test == 'fast':
            df_enrich = compute_geneset_enrichment(args.go_file, args, exclude_ranges, libbgmg_vec, [params], base_params, calc_loglike_diff=True)

            df_go = pd.read_csv(args.go_file_test, sep='\t')
            go_loglike_diff = pd.merge(df_go,
                df_enrich[~df_enrich['GO'].isin(['base', args.go_all_genes_label])][['GO','loglike_diff']].rename(columns={'GO':'GENE'}),
                on='GENE', how='left')[['GO', 'loglike_diff']].groupby('GO').agg(['sum']).reset_index()
            go_loglike_diff.columns = ['_'.join(x).strip('_') for x in go_loglike_diff.columns]
            go_loglike_diff.rename(columns={'loglike_diff_sum':'loglike_diff'}, inplace=True)
            del df_enrich_test['loglike_diff']
            df_enrich_test = pd.merge(df_enrich_test, go_loglike_diff, on='GO', how='left')

        # compute loglike_df accounting for overlapping genes within gene sets
        if args.calc_loglike_diff_go_test in ['fast', 'full']:
            # here args.go_file must be used (not args.go_file_test) because it must be the parameters of the model
            # that count towards df (degrees of freedom)
            go, _ = load_go_file(args.go_file, args.go_extend_bp)
            go = go[~go['GO'].isin(['base', args.go_all_genes_label])]
            tree=intervaltree.IntervalTree.from_tuples(zip(go['FROM'].values, go['TO'].values, go['GO'].values))
            overlaps = []
            for interval_a in tree:
                for interval_b in tree.overlap(interval_a):
                    overlaps.append((interval_a.data, interval_b.data))
            df_overlaps = pd.DataFrame(overlaps, columns=['GENE', 'GENE2'])

            loglike_df = pd.merge(pd.read_csv(args.go_file_test, sep='\t'),
                                df_overlaps,
                                on='GENE')[['GO', 'GENE2']].drop_duplicates().groupby('GO').count().reset_index().rename(columns={'GENE2':'loglike_df'})
            df_enrich_test = pd.merge(df_enrich_test, loglike_df, on='GO', how='left')
            df_enrich_test['loglike_aic'] = -2*df_enrich_test['loglike_df'] + 2*df_enrich_test['loglike_diff']

        df_enrich_test.to_csv(args.out + '.go_test_enrich.csv', sep='\t', index=False)
        libbgmg.log_message("Enrichment results for {} are saved to {}".format(args.go_file_test, args.out + '.go_test_enrich.csv'))

    Path(args.out + '.done').touch()

def execute_fixed_effects_parser(args):
    initialize_libbgmg_logging(args)
    fix_and_validate_args(args)

    bim = pd.concat([pd.read_csv(args.bim_file.replace('@', str(chr_label)), sep='\t', header=None, names='CHR SNP GP BP A1 A2'.split()) for chr_label in args.chr2use])
    del bim['GP']
    bim.reset_index(inplace=True, drop=True)

    causal_variants = pd.read_csv(args.causal_variants, delim_whitespace=True, header=None, names=['SNP', 'BETA'])
    causalbetavec = pd.merge(bim[['SNP']], causal_variants, on='SNP', how='left')[['BETA']].values
    causalbetavec[~np.isfinite(causalbetavec)] = 0
    
    libbgmg = LibBgmg(args.lib)
    if args.loadlib_file:
        libbgmg.load(args.loadlib_file)
    else:
        libbgmg.set_option('use_complete_tag_indices', args.use_complete_tag_indices)
        libbgmg.init(args.bim_file, "", args.chr2use, "", "", args.exclude, args.extract)
        for chr_label in args.chr2use: 
            libbgmg.set_ld_r2_coo_from_file(chr_label, args.ld_file.replace('@', str(chr_label)))
            libbgmg.set_ld_r2_csr(chr_label)
        libbgmg.set_option('diag', 0)
        if args.savelib_file:
            if not args.use_complete_tag_indices: libbgmg.log_message('WARNING: unless you know exactly what you are doing, it is recommended to use --use-complete-tag-indices together with --savelib-file.')
            libbgmg.save(args.savelib_file)
 
    trait_index = 1
    libbgmg.set_causalbetavec(causalbetavec.astype(np.float32).flatten(), trait_index)
    libbgmg.set_nvec(np.ones((libbgmg.num_snp, ), dtype=np.float32), trait_index)
    fixedeffectdelta = libbgmg.get_fixedeffectdelta(trait_index)  # num_tag
    fixedeffectbeta = np.divide(fixedeffectdelta, np.sqrt(libbgmg.get_nvec(trait_index)))

    bim_tag = bim[libbgmg.defvec].copy()
    bim_tag['BETA'] = fixedeffectbeta
    bim_tag['PVAL'] = 1
    bim_tag.to_csv(args.out + '.fixed_effects.csv', sep='\t', index=False)
    libbgmg.log_message('Done')

def execute_split_sumstats_parser(args):
    fix_and_validate_args(args)
    if '@' not in args.out:
        raise(ValueError('--out must contain @ symbol for "mixer.py split_sumstats" call'))
    print(f'reading {args.trait1_file}...')
    df = pd.read_csv(args.trait1_file, sep='\t', dtype=str)
    for chr_label in args.chr2use:
        fout = args.out.replace('@', str(chr_label))
        print(f'writing {fout}...')
        df[df['CHR'] == str(chr_label)].to_csv(fout, sep='\t', index=False)
    print('Done.')

