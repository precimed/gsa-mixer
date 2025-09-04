#!/usr/bin/env python

import argparse
import collections
import json
import numpy as np
import os
import pandas as pd

import scipy.optimize
import scipy.io as sio
import time
import sys

from datetime import datetime
from pathlib import Path
from scipy.sparse import csr_matrix, coo_matrix

from common.libbgmg import LibBgmg, LibBgmgAnnot
from gsa_mixer.utils import AnnotUnivariateParams
import gsa_mixer.cli
from .utils import BivariateParams
from .utils import _params_to_dict
from .utils import _dict_to_params
from .utils import UnivariateParametrization_natural_axis		                     # diffevo, neldermead, uncertainty
from .utils import UnivariateParametrization_constPI_constSIG2BETA                   # inflation
from .utils import UnivariateParametrization_constPI		                         # infinitesimal
from .utils import BivariateParametrization_constUNIVARIATE_constRHOBETA_constPI     # inflation
from .utils import BivariateParametrization_constUNIVARIATE_natural_axis             # diffevo, neldermead
from .utils import BivariateParametrization_pi1_pi2_pi12                             # neldermead-pi
from .utils import BivariateParametrization_constUNIVARIATE_constRG_constRHOZERO     # brute1, brent1
from .utils import BivariateParametrization_constUNIVARIATE_infinitesimal            # infinitesimal

from common.libbgmg import _auxoption_none, _auxoption_tagpdf
from common.libbgmg import _cost_calculator_convolve, _cost_calculator_gaussian, _cost_calculator_sampling, _cost_calculator_smplfast
from common.libbgmg import tagvec_as_snpvec

from .utils import _calculate_univariate_uncertainty
from .utils import _calculate_univariate_uncertainty_funcs
from .utils import _calculate_bivariate_uncertainty_funcs

from .utils_obsolete import UnivariateParams_obsolete

from .utils import calc_qq_plot
from .utils import calc_power_curve
from .utils import calc_bivariate_qq
from .utils import calc_bivariate_pdf

import common.utils_cli

from common.utils import NumpyEncoder

_cost_calculator = {
    'gaussian': _cost_calculator_gaussian, 
    'sampling': _cost_calculator_sampling,
    'convolve': _cost_calculator_convolve,
    'smplfast': _cost_calculator_smplfast,
}

def enhance_optimize_result(r, cost_n, cost_df=None, cost_fast=None):
    # optimize_result - an instance of scipy.optimize.OptimizeResult
    # const_n - number of genetic variants effectively contributing to the cost (sum of weights)
    r['cost_n'] = cost_n
    r['cost_df'] = len(r.x) if (cost_df is None) else cost_df
    r['cost'] = r.fun
    r['BIC'] = np.log(r['cost_n']) * r['cost_df'] + 2 * r['cost']
    r['AIC'] =                   2 * r['cost_df'] + 2 * r['cost']
    if cost_fast is not None: r['cost_fast'] = cost_fast

def parser_add_arguments(parser, func, analysis_type):
    common.utils_cli.parser_add_argument_bim_file(parser)
    common.utils_cli.parser_add_argument_ld_file(parser)
    common.utils_cli.parser_add_argument_chr2use(parser)
    common.utils_cli.parser_add_argument_extract_and_exclude(parser)
    common.utils_cli.parser_add_argument_use_complete_tag_indices(parser)
    common.utils_cli.parser_add_argument_allow_albiguous_snps(parser)
    common.utils_cli.parser_add_argument_savelib_and_loadlib(parser)
    common.utils_cli.parser_add_argument_annot_file(parser)
    common.utils_cli.parser_add_argument_go_file(parser)

    if analysis_type in ['fit1', 'test1', 'fit2', 'test2', 'perf']:
        num_traits = None
        if analysis_type in ['fit1', 'test1']: num_traits = 1
        if analysis_type in ['fit2', 'test2', 'perf']: num_traits = 2

        parser.add_argument("--trait1-file", type=str, default="", help="GWAS summary statistics for the first trait (required argument); for 'save_ldcs_reference' analysis it is recommended to split GWAS summary statistics per chromosome; if this is done then --trait1-file should contain symbol '@', which will be replaced by an actual chromosome label. ")
        if num_traits==2: parser.add_argument("--trait2-file", type=str, default="", help="GWAS summary statistics for the second trait (required argument)")
        parser.add_argument('--z1max', type=float, default=None, help="right-censoring threshold for the first trait (default: %(default)s); recommended setting: '--z1max 9.336' (equivalent to p-value 1e-20)")
        if num_traits==2: parser.add_argument('--z2max', type=float, default=None, help="right-censoring threshold for the second trait (default: %(default)s); recommended setting: '--z2max 9.336'")
        parser.add_argument('--hardprune-subset', type=float, default=0, help="fraction of SNPs to randomly select (default: %(default)s); passing 0 as argument disables the operation; values below 1 are interpreted as fraction, values above 1 are interpreted as number of SNPs to select. Both options consider all SNPs in the reference (not just tag SNPs)")
        parser.add_argument('--hardprune-maf', type=float, default=0.05, help="threshold for minor allele frequency (default: %(default)s); passing 0 as argument disables the operation")
        parser.add_argument('--hardprune-r2', type=float, default=1.0, help="LD r2 threshold for prunning SNPs (default: %(default)s); passing 0 as argument disables the operation")
        parser.add_argument('--randprune-n', type=int, default=0, help="number of random pruning iterations (default: %(default)s)")
        parser.add_argument('--randprune-r2', type=float, default=0.1, help="threshold for random pruning (default: %(default)s)")
        parser.add_argument('--disable-inverse-ld-score-weights', default=False, action="store_true", help="a flag allowing to disable weighting by inverse ld score")

        parser.add_argument('--weights-file', type=str, default=('auto' if (analysis_type in ['fit1', 'fit2']) else 'none'), help="file to load weights (default: %(default)s); 'auto' takes weights file that is associated with --load-params argument (if it is provided); 'none' disables 'auto' feature; otherwise this argument must point to a <out>.weights produced by a previous mixer.py analysis")
        parser.add_argument('--save-weights', default=False, action="store_true", help="a flag allowing to save weights set via --hardprune and --randprune options")

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

    if analysis_type in ['fit1', 'test1', 'fit2', 'test2']:
        parser.add_argument('--nckoef', default=None, type=float, help="(optional) coefficient translating total number of causal variants to the number of causal variants explaining 90%% of trait's heritability; by default this is estimated from the data (up until MiXeR v1.3 this was set to 0.319)")

        if analysis_type in ['fit1', 'test1']:
            parser.add_argument('--cubature-rel-error', type=float, default=1e-5, help="relative error for cubature stop criteria (default: %(default)s); applies to 'convolve' cost calculator")
            parser.add_argument('--cubature-max-evals', type=float, default=1000, help="max evaluations for cubature stop criteria (default: %(default)s); applies to 'convolve' cost calculator. "
                "Bivariate cubature require in the order of 10^4 evaluations and thus is much slower than sampling, therefore it is not exposed via mixer.py command-line interface. ")

        if analysis_type == 'fit1':
            parser.add_argument('--analysis', type=str, default='fit1', help=argparse.SUPPRESS, choices=['fit1', 'fit2', 'test1', 'test2'])
            parser.add_argument('--fit-sequence', type=str, default=['diffevo-fast', 'neldermead'], nargs='+', help=argparse.SUPPRESS, choices=['diffevo', 'diffevo-fast', 'neldermead', 'neldermead-fast', 'infinitesimal'])
        elif analysis_type == 'fit2':
            parser.add_argument('--analysis', type=str, default='fit2', help=argparse.SUPPRESS, choices=['fit1', 'fit2', 'test1', 'test2'])
            parser.add_argument('--fit-sequence', type=str, default=['diffevo-fast', 'neldermead-fast', 'brute1', 'brent1'], nargs='+', help=argparse.SUPPRESS, choices=['diffevo-fast', 'neldermead-fast', 'neldermead-pi', 'neldermead-pi-fast', 'brute1', 'brute1-fast', 'brent1', 'brent1-fast', 'infinitesimal'])
            parser.add_argument('--trait1-params-file', type=str, default=None, help="univariate params for the first trait (required argument); applies to 'fit2' analysis only")
            parser.add_argument('--trait2-params-file', type=str, default=None, help="univariate params for the second trait (required argument); applies to 'fit2' analysis only")
        elif analysis_type == 'test1':
            parser.add_argument('--analysis', type=str, default='test1', help=argparse.SUPPRESS, choices=['fit1', 'fit2', 'test1', 'test2'])
            parser.add_argument('--fit-sequence', type=str, default=['inflation'], nargs='+', help=argparse.SUPPRESS, choices=['inflation'])
        elif analysis_type == 'test2':
            parser.add_argument('--analysis', type=str, default='test2', help=argparse.SUPPRESS, choices=['fit1', 'fit2', 'test1', 'test2'])
            parser.add_argument('--fit-sequence', type=str, default=['inflation'], nargs='+', help=argparse.SUPPRESS, choices=['inflation'])
        else:
            raise ValueError('internal error: invalid combination of do_fit and num_traits')

        if analysis_type in ['fit1', 'fit2', 'test1', 'test2']:
            parser.add_argument('--load-params-file', type=str, default=None, help="params of the fitted model (required argument for test1 and test2, otherwise optional);")

        # all arguments below marked with argparse.SUPRESS option are internal and not recommended for a general use.
        # Valid options for --fit-sequence (remember to combine them in a sequence that makes sence):
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

    parser.set_defaults(func=func)

def load_params_file(fname_or_params, libbgmg, args):
    if isinstance(fname_or_params, AnnotUnivariateParams) or isinstance(fname_or_params, BivariateParams):
        return fname_or_params  # used for unit tests
    if isinstance(fname_or_params, str):
        data = json.loads(open(fname_or_params).read())
        return _dict_to_params(data['params'], libbgmg, args)
    else:
        return _dict_to_params(fname_or_params, libbgmg, args)

def apply_univariate_fit_sequence(args, libbgmg, fit_sequence, init_params=None, trait=1):
    params=init_params if (init_params is not None) else AnnotUnivariateParams(pi=1.0, libbgmg=libbgmg, s=0, l=0, sig2_zeroL=0)
    optimize_result_sequence=[]

    if (not init_params) and args.load_params_file:
        libbgmg.log_message("Loading {}...".format(args.load_params_file))
        params = load_params_file(args.load_params_file, libbgmg, args)
        libbgmg.log_message("Done, {}".format(params))

    for fit_type in fit_sequence:
        optimize_result = None
        libbgmg.log_message("fit_type=={}...".format(fit_type))

        if (fit_type == 'diffevo') or fit_type == ('diffevo-fast'):
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
    params = None; optimize_result_sequence = []

    if args.load_params_file:
        libbgmg.log_message("Loading initial bivariate params {}...".format(args.load_params_file))
        params = load_params_file(args.load_params_file, libbgmg, args)
        libbgmg.log_message("Done, {}".format(params))
        params1 = params._params1
        params2 = params._params2
    else:
        libbgmg.log_message("Loading univariate constrains for bivariate analysis...")
        params1 = load_params_file(args.trait1_params_file, libbgmg, args)
        params2 = load_params_file(args.trait2_params_file, libbgmg, args)
        libbgmg.log_message("trait1: {}".format(params1))
        libbgmg.log_message("trait2: {}".format(params2))

    for fit_type in fit_sequence:
        optimize_result = None
        libbgmg.log_message("fit_type=={}...".format(fit_type))
        if fit_type == 'inflation':
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

        elif (fit_type == 'neldermead-pi') or (fit_type == 'neldermead-pi-fast'):
            if params == None: raise(RuntimeError('params == None, unable to proceed with "{}" fit'.format(fit_type)))
            libbgmg.set_option('cost_calculator', _cost_calculator_sampling if (fit_type == 'neldermead-pi') else _cost_calculator_gaussian)
            cost_prior_to_optimization = params.cost(libbgmg)
            parametrization = BivariateParametrization_pi1_pi2_pi12(params=params, lib=libbgmg)
            optimize_result = scipy.optimize.minimize(lambda x: parametrization.calc_cost(x), parametrization.params_to_vec(params),
                method='Nelder-Mead', options={'maxiter':1200, 'fatol':1e-7, 'xatol':1e-4, 'adaptive':True})
            params = parametrization.vec_to_params(optimize_result.x)
            optimize_result['cost-prior-to-optimization'] =cost_prior_to_optimization

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
        if optimize_result:
            enhance_optimize_result(optimize_result, cost_df=9, cost_n=np.sum(libbgmg.weights), cost_fast=params.cost(libbgmg))
            optimize_result['params']=_params_to_dict(params)   # params after optimization
            optimize_result_sequence.append((fit_type, optimize_result))
        libbgmg.log_message("fit_type=={} done ({}, {})".format(fit_type, params, optimize_result))

    if params == None: raise(RuntimeError('Empty --fit-sequence'))
    return (params, params1, params2, optimize_result_sequence)

def execute_perf_parser(args):
    common.utils_cli.fix_and_validate_args(args)
    libbgmg = gsa_mixer.cli.initialize_mixer_plugin(args, 0, args.trait1_file, args.trait2_file, args.randprune_n, args.loadlib_file, args.savelib_file, args.chr2use, '@')
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
    common.utils_cli.fix_and_validate_args(args)

    libbgmg = gsa_mixer.cli.initialize_mixer_plugin(args, 0, args.trait1_file, "", args.randprune_n, args.loadlib_file, args.savelib_file, args.chr2use, '@')

    if not args.trait1_file:
        libbgmg.log_message('--trait1-file not provided, exit')
        return

    df_weights = pd.DataFrame({'chrnumvec':libbgmg.chrnumvec, 'weights': tagvec_as_snpvec(libbgmg.weights, libbgmg.defvec)})
    results = init_results_struct(args, [libbgmg])
    results['analysis'] = 'univariate'

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
    common.utils_cli.fix_and_validate_args(args)

    libbgmg = gsa_mixer.cli.initialize_mixer_plugin(args, 0, args.trait1_file, args.trait2_file, args.randprune_n, args.loadlib_file, args.savelib_file, args.chr2use, '@')

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

def execute_save_matlab_reference_parser(args):
    common.utils_cli.fix_and_validate_args(args)

    libbgmg_vec = []
    annot_snps = 0; gene_snps = 0
    for chr_label in args.chr2use:
        libbgmg = gsa_mixer.cli.initialize_mixer_plugin(
            args,
            context_id=chr_label,
            trait1_file="",
            trait2_file="",
            randprune_n=None, 
            loadlib_file=args.loadlib_file.replace('@', str(chr_label)) if args.loadlib_file else None, 
            savelib_file=args.savelib_file.replace('@', str(chr_label)) if args.savelib_file else None, 
            chr_labels=[chr_label], 
            lib_file_suffix=str(chr_label))
        annonames = libbgmg.annonames
        genenames = libbgmg.genenames
        libbgmg_vec.append(libbgmg)
        annot_snps = annot_snps + np.sum(libbgmg.annomat, 0).flatten()
        gene_snps = gene_snps + np.sum(libbgmg.genemat, 0).flatten()

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

    Path(args.out + '.done').touch()

def execute_fixed_effects_parser(args):
    common.utils_cli.fix_and_validate_args(args)

    bim = pd.concat([pd.read_csv(args.bim_file.replace('@', str(chr_label)), sep='\t', header=None, names='CHR SNP GP BP A1 A2'.split()) for chr_label in args.chr2use])
    del bim['GP']
    bim.reset_index(inplace=True, drop=True)

    causal_variants = pd.read_csv(args.causal_variants, sep=r'\s+', header=None, names=['SNP', 'BETA'])
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
