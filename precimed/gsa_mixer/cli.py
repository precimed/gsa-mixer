#!/usr/bin/env python

import argparse
import collections
import json
import numpy as np
import os
import pandas as pd
import random
import scipy.io as sio
import sys
import intervaltree
import scipy

from datetime import datetime
from pathlib import Path
from scipy.sparse import csr_matrix, coo_matrix

import common.utils_cli
from common.libbgmg import LibBgmg, LibBgmgAnnot
from common.libbgmg import tagvec_as_snpvec

from .utils import AnnotUnivariateParams
from .utils import _params_to_dict
from .utils import _dict_to_params
from .utils import AnnotUnivariateParametrization

from common.utils import NumpyEncoder
from common.libbgmg import _cost_calculator_sampling

CHR_OFFSET = np.int64(1e10)

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
            common.utils_cli.check_input_file(args, 'weights_file')

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

def parser_plsa_add_arguments(parser, func):
    common.utils_cli.parser_add_argument_bim_file(parser)
    common.utils_cli.parser_add_argument_ld_file(parser)
    common.utils_cli.parser_add_argument_chr2use(parser)
    common.utils_cli.parser_add_argument_extract_and_exclude(parser)
    common.utils_cli.parser_add_argument_use_complete_tag_indices(parser)
    common.utils_cli.parser_add_argument_allow_albiguous_snps(parser)
    common.utils_cli.parser_add_argument_savelib_and_loadlib(parser)

    parser.add_argument("--trait1-file", type=str, default="", help="GWAS summary statistics for the first trait (required argument); for 'plsa' analysis it is recommended to split GWAS summary statistics per chromosome; if this is done then --trait1-file should contain symbol '@', which will be replaced by an actual chromosome label. ")
    parser.add_argument('--z1max', type=float, default=None, help="right-censoring threshold for the first trait (default: %(default)s); recommended setting: '--z1max 9.336' (equivalent to p-value 1e-20)")
    parser.add_argument('--hardprune-subset', type=float, default=0, help="fraction of SNPs to randomly select (default: %(default)s); passing 0 as argument disables the operation; values below 1 are interpreted as fraction, values above 1 are interpreted as number of SNPs to select. Both options consider all SNPs in the reference (not just tag SNPs)")
    parser.add_argument('--hardprune-maf', type=float, default=0.05, help="threshold for minor allele frequency (default: %(default)s); passing 0 as argument disables the operation")
    parser.add_argument('--hardprune-r2', type=float, default=1.0, help="LD r2 threshold for prunning SNPs (default: %(default)s); passing 0 as argument disables the operation")
    parser.add_argument('--randprune-n', type=int, default=0, help="number of random pruning iterations (default: %(default)s)")
    parser.add_argument('--randprune-r2', type=float, default=0.1, help="threshold for random pruning (default: %(default)s)")
    parser.add_argument('--disable-inverse-ld-score-weights', default=False, action="store_true", help="a flag allowing to disable weighting by inverse ld score")

    parser.add_argument('--weights-file', type=str, default='auto', help="file to load weights (default: %(default)s); 'auto' takes weights file that is associated with --load-params argument (if it is provided); 'none' disables 'auto' feature; otherwise this argument must point to a <out>.weights produced by a previous mixer.py analysis")
    parser.add_argument('--save-weights', default=False, action="store_true", help="a flag allowing to save weights set via --hardprune and --randprune options")

    parser.add_argument('--seed', type=np.uint32, default=None, help="random seed (default: %(default)s).")
    parser.add_argument('--r2min', type=float, default=0.0, help=argparse.SUPPRESS) # help="r2 values below this threshold will contribute via infinitesimal model")
    parser.add_argument('--threads', type=int, default=[None], nargs='+', help="specify how many concurrent threads to use for computations; (default: total number of CPU cores available)")

    parser.add_argument('--kmax', type=int, default=[20000], nargs='+', help="number of sampling iterations for log-likelihod and posterior delta (default: 20000)")

    common.utils_cli.parser_add_argument_annot_file(parser, annot_file_test=True)
    common.utils_cli.parser_add_argument_go_file(parser, go_file_test=True)
    parser.add_argument("--calc-loglike-diff-go-test", default=None, choices=['fast', 'full'], help="algorithm for computation of the loglike_diff column in <out>.go_test_enrich.csv")
    parser.add_argument("--go-all-genes-label", type=str, default='coding_genes', help="reference gene-set to calibrate fold enrichment, e.g. allowing to compute enrichment w.r.t. the set of all coding genes (default: %(default)s)")
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
    parser.add_argument('--adam-disable', default=False, action="store_true", help="a flag allowing to disable optimization; require '--load-params-file <out-of-a-previous-run>.json'")
    parser.add_argument('--adam-separate-by-chromosome', default=False, action="store_true", help="a flag allowing to run optimization separately on each chromosome; this is only recommended with '--fit gene --constrain-base-gene-category --adam-beta2 0.9' configuration")
    
    parser.add_argument('--load-params-file', type=str, default=None, help="(optional) params of the fitted model to load; expected to be from a 'mixer.py plsa' run; if specified, parameters loaded from --load-params-file take priority over --s-value, --l-value, --pi-value, --sig2-zeroA, --sig2-zeroL values; note that --go-all-genes-label is used to populate initial sig2_gene values for all genes except base")
    parser.add_argument('--load-baseline-params-file', type=str, default=None, help="(optional) params of the baseline model to load and use for fold enrichment estimation; if --load-baseline-params-file is not provided, --load-params-file will be used;")

    parser.add_argument("--calc-detailed-trajectory", default=False, action="store_true", help=argparse.SUPPRESS) # "Calcualte detailed trajectory (one data point for each chromosome and --adam-epoch). Make cause excessively large .json files with --go-annot.")
    parser.add_argument("--calc-trajectory", default=False, action="store_true", help=argparse.SUPPRESS)
    parser.add_argument('--se-samples', default=100, type=int, help="number of samples in standard error computation (using log-likelihood hessian / observed Fisher's information; applies only to --go-file-test")
    parser.add_argument("--gsa-base", default=False, action="store_true", help="set default parameters to values recommended for the baseline GSA model")
    parser.add_argument("--gsa-full", default=False, action="store_true", help="set default parameters to values recommended for the full GSA model")

    parser.set_defaults(func=func)

def parser_split_sumstats_add_arguments(func, parser):
    parser.add_argument("--trait1-file", type=str, default="", help="GWAS summary statistics for the first trait (required argument);")
    common.utils_cli.parser_add_argument_chr2use(parser)
    parser.set_defaults(func=func)

def parser_ld_add_arguments(func, parser):
    parser.add_argument("--bfile", type=str, default=None, help="path to plink bfile (required argument)")
    parser.add_argument('--r2min', type=float, default=0.05, help="r2 values above this threshold will be stored in sparse LD format (default: %(default)s)")
    parser.add_argument('--ldscore-r2min', type=float, default=0.001, help="r2 values above this threshold (and below --r2min) will be stored as LD scores that contribute to the cost function via an infinitesimal model (default: %(default)s)")
    parser.add_argument('--ld-window-kb', type=float, default=0, help="limit window similar to --ld-window-kb in 'plink r2'; 0 will disable this constraint (default: %(default)s); either ld-window-kb or ld-window argument must be provided")
    parser.add_argument('--ld-window', type=int, default=0, help="limit window similar to --ld-window in 'plink r2'; 0 will disable this constraint (default: %(default)s); either ld-window-kb or ld-window argument must be provided")
    parser.set_defaults(func=func)

def load_params_file(fname_or_params, libbgmg, args):
    if isinstance(fname_or_params, AnnotUnivariateParams):
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

def execute_ld_parser(args):
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
    exclude_ranges = ' '.join([f'{range.chr}:{int(range.from_bp)}-{int(range.to_bp)}' for range in common.utils_cli.make_ranges(args.exclude_ranges)]) if ('exclude_ranges' in args) else ""

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

        for opt, val in common.utils_cli.convert_args_to_libbgmg_options(args):
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

        for opt, val in common.utils_cli.convert_args_to_libbgmg_options(args):
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

def find_per_snp_information_univariate(libbgmg_vec, params):
    # Save a file with the following information
    df_snps = []
    for libbgmg in libbgmg_vec:
        params._libbgmg = libbgmg
        params._snp_defvec = None
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
    libbgmg = LibBgmg(args.lib)
    common.utils_cli.fix_and_validate_args(args)
    if (args.load_params_file is None) and (args.adam_disable):
        raise ValueError('--load-params-file is required for --adam-disable')

    if args.go_all_genes_label == 'base': raise ValueError('"--go-all-genes-label base" is not allowed')

    if args.gsa_base:
        args.fit = ['sig2-zeroA', 'l', 'annot', 'gene']
        args.l_value = -0.25
        args.h2_init_calibration = 0.1
        args.se_samples = 0
        args.save_weights = True
        libbgmg.log_message(f'apply --gsa-base settings (--fit {" ".join(args.fit)} --l-value {args.l_value} --h2-init-calibration {args.h2_init_calibration}) --save-weights ')

    if args.gsa_full:
        args.fit = ['gene']
        args.constrain_base_gene_category = True
        args.nullify_go_all_genes_sig2_gene = True
        args.calc_loglike_diff_go_test = 'fast'
        args.adam_separate_by_chromosome = True
        args.adam_beta1 = 0.1
        args.adam_beta2 = 0.8
        libbgmg.log_message(f'apply --gsa-full settings (--fit {" ".join(args.fit)}  --constrain-base-gene-category --nullify-go-all-genes-sig2-gene --calc-loglike-diff-go-test fast --adam-separate-by-chromosome --adam-beta1 {args.adam_beta1} --adam-beta2 {args.adam_beta2})')

    if args.annot_file_test or args.go_file_test:
        if (not args.load_params_file) and (not args.load_baseline_params_file):
            raise ValueError('--annot-file-test/--go-file-test require either --load-params-file or --load-baseline-params-file to be provided')
        if (args.load_baseline_params_file is None) and (args.load_params_file is not None):
            args.load_baseline_params_file = args.load_params_file
        common.utils_cli.check_input_file(args, 'load_baseline_params_file')

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

    if args.gsa_base:
        libbgmg.log_message('Finding per-SNP information...')
        df=find_per_snp_information_univariate(libbgmg_vec, params)
        fname = args.out + '.snps.csv'
        df.to_csv(fname, sep='\t', index=False)
        libbgmg.log_message('Saved {} file covering {} SNPs'.format(fname, len(df)))

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
    exclude_ranges = common.utils_cli.make_ranges(args.exclude_ranges)
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

def execute_split_sumstats_parser(args):
    common.utils_cli.fix_and_validate_args(args)
    if '@' not in args.out:
        raise(ValueError('--out must contain @ symbol for "mixer.py split_sumstats" call'))
    print(f'reading {args.trait1_file}...')
    df = pd.read_csv(args.trait1_file, sep='\t', dtype=str)
    for chr_label in args.chr2use:
        fout = args.out.replace('@', str(chr_label))
        print(f'writing {fout}...')
        df[df['CHR'] == str(chr_label)].to_csv(fout, sep='\t', index=False)
    print('Done.')

