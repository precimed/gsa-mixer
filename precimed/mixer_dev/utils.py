# Utility classes for univariate and bivariate fit
# Contains
# _log_exp_converter, _logit_logistic_converter, _arctanh_tanh_converter - converters to map bounded parameters into -inf, +inf range
# UnivariateParams, BivariateParams - represent parameters, with basic functionality like "calculate cost"
# Several univariate and bivariate parametrizations - suitable for fitting and (some specific parametrization) for uncertainty calculation
# _calculate_univariate_uncertainty - estimates confidence intervals for parameters and their aggregates
 
import numpy as np
import numdifftools as nd
import scipy.stats
from scipy.interpolate import interp1d
from common.utils import _arctanh_tanh_converter, _log_exp_converter, _logit_logistic_converter
from common.libbgmg import _auxoption_none, _auxoption_tagpdf

import gsa_mixer.utils

def find_rg_sig2_factor(params1, params2):
    lib = params1._libbgmg
    hetvec = np.float32(2.0) * lib._mafvec * (1-lib._mafvec)
    het_sig2_trait1 = np.multiply(hetvec, params1.find_sig2_mat().flatten())
    het_sig2_trait2 = np.multiply(hetvec, params2.find_sig2_mat().flatten())
    het_cov = np.sqrt(np.multiply(het_sig2_trait1, het_sig2_trait2))
    return np.sum(het_cov) / np.sqrt(np.sum(het_sig2_trait1) * np.sum(het_sig2_trait2))

class BivariateParams(object):
    def __init__(self, params1=None, params2=None, rho_beta=None, rho_zero=None, pi12=None):
        self._params1 = params1
        self._params2 = params2

        self._pi12 = pi12
        self._rho_beta = rho_beta
        self._rho_zero = rho_zero

        self._validate()

    def _rg_sig2_factor(self):
        return find_rg_sig2_factor(self._params1, self._params2)

    def _rg(self):
        return self._rg_sig2_factor() * self._rho_beta * self._pi12 / np.sqrt(self._params1.pi * self._params2.pi)

    def find_nckoef(self):
        # NCKoef is a coefficient between 0 and 1 that represents the proportion of causal variants explaining 90% of trait's heritability.
        # In MiXeR v1.3 this was hard-coded to 0.319, a value specific to the basic MiXeR model (1-pi1)*N(0,0) + pi1*N(0, sigma2_beta),
        # and estimated from 1kG EUR reference panel.
        # In MiXeR v2.0 and onwards this coefficient depends on other parameters of the AnnotUnivariateParams models,
        # and its also evaluated based on the reference panel being used.
        NCKoef1 = self._params1.find_nckoef()
        NCKoef2 = self._params2.find_nckoef()
        return np.sqrt(NCKoef1 * NCKoef2)

    def _validate(self):
        assert(self._params1._sig2_zeroL == 0)
        assert(self._params2._sig2_zeroL == 0)

        assert np.isscalar(self._pi12)
        assert np.isscalar(self._rho_zero)
        assert np.isscalar(self._rho_beta)
        
        assert (self._pi12 >= 0) and (self._pi12 <= self._params1.pi) and (self._pi12 <= self._params2.pi)
        assert np.greater_equal([self._rho_zero, self._rho_beta], -1).all()
        assert np.less_equal([self._rho_zero, self._rho_beta], 1).all()

    def __str__(self):
        description = []
        for attr_name in '_pi12', '_rho_beta', '_rho_zero':
            try:
                attr_value = getattr(self, attr_name)
                description.append('{}: {}'.format(attr_name, attr_value))
            except RuntimeError:
                pass
        description.append('rg_sig2_factor: {}'.format(self._rg_sig2_factor()))
        description.append('rg: {}'.format(self._rg()))
        return 'BivariateParams({})'.format(', '.join(description))
    __repr__ = __str__
    
    def find_pi_mat(self, num_snp):
        pi = [self._params1.pi - self._pi12, self._params2.pi - self._pi12, self._pi12]
        return np.matmul(np.array(pi, dtype=np.float32).reshape(-1, 1), np.ones(shape=(1, num_snp), dtype=np.float32))

    def find_sig2_mat(self):
        return np.concatenate((np.transpose(self._params1.find_sig2_mat()), np.transpose(self._params2.find_sig2_mat())), axis=0)

    def find_rho_vec(self, num_snp):
        return self._rho_beta * np.ones(shape=(num_snp, 1), dtype=np.float32)

    def cost(self, lib):
        num_snp = lib.num_snp
        value = lib.calc_unified_bivariate_cost(self.find_pi_mat(num_snp), self.find_sig2_mat(), self.find_rho_vec(num_snp),
                                                sig2_zeroA=[self._params1._sig2_zeroA, self._params2._sig2_zeroA],
                                                sig2_zeroC=[1, 1], sig2_zeroL=[0, 0],
                                                rho_zeroA=self._rho_zero,
                                                rho_zeroL=0)
        return value if np.isfinite(value) else 1e100

    def aux(self, lib):
        num_snp = lib.num_snp
        return lib.calc_unified_bivariate_aux(self.find_pi_mat(num_snp), self.find_sig2_mat(), self.find_rho_vec(num_snp),
                                              sig2_zeroA=[self._params1._sig2_zeroA, self._params2._sig2_zeroA],
                                              sig2_zeroC=[1, 1], sig2_zeroL=[0, 0],
                                              rho_zeroA=self._rho_zero,
                                              rho_zeroL=0)

    def tag_pdf(self, lib):
        lib.set_option('aux_option', _auxoption_tagpdf)
        retval = self.aux(lib)
        lib.set_option('aux_option', _auxoption_none)
        return np.array(retval.flat[:lib.num_tag])

    def pdf(self, lib, zgrid):
        num_snp = lib.num_snp
        [zgrid1, zgrid2] = np.meshgrid(zgrid, zgrid)
        zgrid1=zgrid1[zgrid>=0, :]; zgrid2=zgrid2[zgrid>=0, :]

        pdf = lib.calc_unified_bivariate_pdf(self.find_pi_mat(num_snp), self.find_sig2_mat(), self.find_rho_vec(num_snp),
                                             sig2_zeroA=[self._params1._sig2_zeroA, self._params2._sig2_zeroA],
                                             sig2_zeroC=[1, 1], sig2_zeroL=[0, 0],
                                             rho_zeroA=self._rho_zero,
                                             rho_zeroL=0,
                                             zvec1=zgrid1.flatten(), zvec2=zgrid2.flatten())

        pdf = pdf.reshape(zgrid1.shape)
        pdf = np.concatenate((np.fliplr(np.flipud(pdf[1:, :])), pdf))
        return pdf

    def delta_posterior(self, lib):
        num_snp = lib.num_snp
        return lib.calc_unified_bivariate_delta_posterior(
            self.find_pi_mat(num_snp), self.find_sig2_mat(), self.find_rho_vec(num_snp),
            sig2_zeroA=[self._params1._sig2_zeroA, self._params2._sig2_zeroA],
            sig2_zeroC=[1, 1], sig2_zeroL=[0, 0],
            rho_zeroA=self._rho_zero,
            rho_zeroL=0)

# Parametrization that constrains polygenicity
class UnivariateParametrization_constPI(object):
    def __init__(self, params, lib, trait):
        self._params = params.copy()
        self._lib = lib
        self._trait = trait

    def vec_to_params(self, vec):
        self._params._sig2_beta = _log_exp_converter(vec[0], invflag=True)
        self._params._sig2_zeroA = _log_exp_converter(vec[1], invflag=True)
        return self._params

    def params_to_vec(self, params):
        return [_log_exp_converter(params.sig2_beta, invflag=False),
                _log_exp_converter(params.sig2_zeroA, invflag=False)] 

    def calc_cost(self, vec):
        return self.vec_to_params(vec).cost(self._lib, self._trait)

# Parametrization that constraints pi and sig2_beta parameters
class UnivariateParametrization_constPI_constSIG2BETA(object):
    def __init__(self, params, lib, trait):
        self._params = params.copy()
        self._lib = lib    
        self._trait = trait

    def vec_to_params(self, vec):
        self._params._sig2_zeroA = _log_exp_converter(vec[0], invflag=True)
        return self._params
    
    def params_to_vec(self, params):
        return [_log_exp_converter(params.sig2_zeroA, invflag=False)]

    def calc_cost(self, vec):
        return self.vec_to_params(vec).cost(self._lib, self._trait)

# Parametrization that constraints pi and sig2_zero parameters
class UnivariateParametrization_constPI_constSIG2ZERO(object):
    def __init__(self, params, lib, trait):
        self._params = params.copy()
        self._lib = lib
        self._trait = trait

    def vec_to_params(self, vec):
        self._params._sig2_beta = _log_exp_converter(vec[0], invflag=True)
        return self._params

    def params_to_vec(self, params):
        return [_log_exp_converter(params.sig2_beta, invflag=False)]

    def calc_cost(self, vec):
        return self.vec_to_params(vec).cost(self._lib, self._trait)

# Standard univariate parametrization
#   x1 = log(sig2_zeroA)
#   x2 = log(sig2_beta)
#   x3 = logit(pi)
class UnivariateParametrization_natural_axis(object):
    def __init__(self, params, lib, trait):
        self._params = params.copy()
        self._lib = lib
        self._trait = trait

    def params_to_vec(self, params):
        return [_log_exp_converter(params._sig2_zeroA, invflag=False),
                _log_exp_converter(params._sig2_beta, invflag=False),
                _logit_logistic_converter(params.pi, invflag=False)]
        
    def vec_to_params(self, vec):
        self._params._sig2_zeroA=_log_exp_converter(vec[0], invflag=True)
        self._params._sig2_beta=_log_exp_converter(vec[1], invflag=True)
        self._params._pi=_logit_logistic_converter(vec[2], invflag=True)
        return self._params
    
    def calc_cost(self, vec):
        return self.vec_to_params(vec).cost(self._lib, self._trait)

class BivariateParametrization_pi1_pi2_pi12(object):
    def __init__(self, params, lib):
        self._lib = lib
        self._params1 = params._params1.copy()
        self._params2 = params._params2.copy()

        self._const_rho_beta = params._rho_beta
        self._const_rho_zero = params._rho_zero

        # Compensate for "sig2_beta_effective = self.sig2_beta / self.pi" from AnnotUnivariateParams.find_sig2_mat
        self._const_effective_sig2_beta1 = self._params1._sig2_beta / self._params1._pi
        self._const_effective_sig2_beta2 = self._params2._sig2_beta / self._params2._pi

    def params_to_vec(self, params):
        pi12 = params._pi12
        pi1 = params._params1._pi - pi12
        pi2 = params._params2._pi - pi12
        return [_logit_logistic_converter(pi1, invflag=False),
                _logit_logistic_converter(pi2, invflag=False),
                _logit_logistic_converter(pi12, invflag=False)]

    def vec_to_params(self, vec):
        pi1 = _logit_logistic_converter(vec[0], invflag=True)
        pi2 = _logit_logistic_converter(vec[1], invflag=True)
        pi12 =_logit_logistic_converter(vec[2], invflag=True)

        # fix univariate pi & sig2_beta
        self._params1._pi = pi1 + pi12
        self._params2._pi = pi2 + pi12
        self._params1._sig2_beta = self._const_effective_sig2_beta1 * (pi1 + pi12)
        self._params2._sig2_beta = self._const_effective_sig2_beta2 * (pi2 + pi12)

        return BivariateParams(params1=self._params1, params2=self._params2,
                               rho_beta=self._const_rho_beta, rho_zero=self._const_rho_zero, pi12=pi12)

    def calc_cost(self, vec):
        return self.vec_to_params(vec).cost(self._lib)

# BGMG_cpp_fit_bivariate_fast (fits pi12 and rho_beta constrained on rg and all univariate params)
# The difference from the above function (....boundedPI) is that here we map params into [-inf,+inf],
# while ...boundedPI is parametrized directly with PI.
class BivariateParametrization_constUNIVARIATE_constRG_constRHOZERO(object):
    def __init__(self, const_params1, const_params2, const_rg, const_rho_zero, lib):
        self._const_params1 = const_params1
        self._const_params2 = const_params2
        self._const_rg = const_rg
        self._const_rg_sig2_factor = find_rg_sig2_factor(self._const_params1, self._const_params2)
        self._const_rho_zero = const_rho_zero
        self._lib = lib

    @property
    def _max_pi12(self):
        return min([self._const_params1.pi, self._const_params2.pi])

    @property
    def _min_pi12(self):
        return abs(self._const_rg) * np.sqrt(self._const_params1.pi * self._const_params2.pi) / self._const_rg_sig2_factor

    def params_to_vec(self, params):
        assert(params._pi12 >= self._min_pi12)
        assert(params._pi12 <= self._max_pi12)
        return [_logit_logistic_converter((params._pi12 - self._min_pi12) / (self._max_pi12 - self._min_pi12), invflag=False)]

    def vec_to_params(self, vec):
        if not hasattr(vec, "__len__"): vec = [vec]
        pi12 = self._min_pi12 + _logit_logistic_converter(vec[0], invflag=True) * (self._max_pi12 - self._min_pi12) 
        rho_beta = self._const_rg * np.sqrt(self._const_params1.pi * self._const_params2.pi) / (pi12 * self._const_rg_sig2_factor)
        assert abs(rho_beta) <= 1

        return BivariateParams(params1=self._const_params1, params2=self._const_params2,
                               rho_beta=rho_beta, rho_zero=self._const_rho_zero, pi12=pi12)

    def calc_cost(self, vec):
        return self.vec_to_params(vec).cost(self._lib)

class BivariateParametrization_constUNIVARIATE_natural_axis(object):
    def __init__(self, const_params1, const_params2, lib):
        self._const_params1 = const_params1
        self._const_params2 = const_params2
        self._lib = lib

    def params_to_vec(self, params):
        max_pi12 = min(self._const_params1.pi, self._const_params2.pi)
        return [_arctanh_tanh_converter(params._rho_beta, invflag=False),
                _arctanh_tanh_converter(params._rho_zero, invflag=False),
                _logit_logistic_converter(params._pi12 / max_pi12, invflag=False)]

    def vec_to_params(self, vec):
        max_pi12 = min(self._const_params1.pi, self._const_params2.pi)
        rho_beta = _arctanh_tanh_converter(vec[0], invflag=True)
        rho_zero = _arctanh_tanh_converter(vec[1], invflag=True)
        pi12 = max_pi12 * _logit_logistic_converter(vec[2], invflag=True)
        return BivariateParams(params1=self._const_params1, params2=self._const_params2,
                               rho_beta=rho_beta, rho_zero=rho_zero, pi12=pi12)

    def calc_cost(self, vec):
        return self.vec_to_params(vec).cost(self._lib)

class BivariateParametrization_constUNIVARIATE_infinitesimal(object):
    def __init__(self, const_params1, const_params2, lib):
        assert(const_params1._pi == 1.0)
        assert(const_params2._pi == 1.0)
        self._const_params1 = const_params1
        self._const_params2 = const_params2
        self._lib = lib

    def params_to_vec(self, params):
        return [_arctanh_tanh_converter(params._rho_beta, invflag=False),
                _arctanh_tanh_converter(params._rho_zero, invflag=False)]

    def vec_to_params(self, vec):
        rho_beta = _arctanh_tanh_converter(vec[0], invflag=True)
        rho_zero = _arctanh_tanh_converter(vec[1], invflag=True)
        return BivariateParams(params1=self._const_params1, params2=self._const_params2,
                               rho_beta=rho_beta, rho_zero=rho_zero, pi12=1.0)

    def calc_cost(self, vec):
        return self.vec_to_params(vec).cost(self._lib)

class BivariateParametrization_constUNIVARIATE_constRHOBETA_constPI(object):
    def __init__(self, const_params1, const_params2, const_pi12, const_rho_beta, lib):
        assert (const_pi12 >= 0) and (const_pi12 <= min(const_params1.pi, const_params2.pi))
        self._const_pi12 = const_pi12
        self._const_rho_beta = const_rho_beta
        self._const_params1 = const_params1
        self._const_params2 = const_params2
        self._lib = lib

    def vec_to_params(self, vec):
        return BivariateParams(params1=self._const_params1, params2=self._const_params2,
                               rho_beta=self._const_rho_beta, rho_zero=_arctanh_tanh_converter(vec[0], invflag=True), pi12=self._const_pi12)

    def params_to_vec(self, params):
        assert abs(params._rho_zero) <= 1
        return [_arctanh_tanh_converter(params._rho_zero, invflag=False)]

    def calc_cost(self, vec):
        return self.vec_to_params(vec).cost(self._lib)
    
    def fit(self, optimizer):
        result = optimizer(self._calc_cost, self._init_vec)
        return self._vec_to_params(result.x), result

def _hessian_robust(hessian, hessdiag):
    # for noisy functions hessian might be badly estimated
    # if we detect a problem with hessian fall back to hess diagonal
    hessdiag[hessdiag < 0] = 1e15
    hessdiag = np.diag(hessdiag)
    if not np.isfinite(hessian).all(): return hessdiag
    try:
        hessinv = np.linalg.inv(hessian)
    except np.linalg.LinAlgError as err:
        return hessdiag
    if not np.isfinite(hessinv).all(): return hessdiag
    if np.less_equal(np.linalg.eigvals(hessinv), 0).any(): return hessdiag
    return hessian

def _calculate_univariate_uncertainty_funcs(alpha, NCKoef):
    funcs = [('pi', lambda x: x.pi),
             ('nc', lambda x: x.pi * x._libbgmg.num_snp),
             ('nc@p9', lambda x: NCKoef * x.pi * x._libbgmg.num_snp),
             ('sig2_beta', lambda x: x.find_mean_sig2_beta()),
             ('sig2_zeroA', lambda x: x._sig2_zeroA),
             ('s', lambda x: x._s),
             ('l', lambda x: x._l),
             ('h2', lambda x: x.find_h2())]
    stats = [('mean', lambda x: np.mean(x)),
             ('median', lambda x: np.median(x)),
             ('std', lambda x: np.std(x)),
             ('lower', lambda x: np.percentile(x, 100.0 * (  alpha/2))),
             ('upper', lambda x: np.percentile(x, 100.0 * (1-alpha/2)))]
    return funcs, stats

def _calculate_univariate_uncertainty(params, parametrization, alpha, num_samples, NCKoef):
    init_vec = parametrization.params_to_vec(params)
    funcs, stats = _calculate_univariate_uncertainty_funcs(alpha, NCKoef)
    hessian = _hessian_robust(nd.Hessian(parametrization.calc_cost)(init_vec), 
                              nd.Hessdiag(parametrization.calc_cost)(init_vec))
    x_sample = np.random.multivariate_normal(init_vec, np.linalg.inv(hessian), num_samples)
    sample = [parametrization.vec_to_params(x) for x in x_sample]
    result = {}
    for func_name, func in funcs:
        result[func_name] = {'point_estimate': func(parametrization.vec_to_params(init_vec))}
        param_vector = [func(s) for s in sample]
        for stat_name, stat in stats:
            result[func_name][stat_name] = stat(param_vector)
    return result, sample

def _calculate_bivariate_uncertainty_funcs(alpha, NCKoef):
    funcs = [('sig2_zero_T1', lambda x: x._params1._sig2_zeroA),
             ('sig2_zero_T2', lambda x: x._params2._sig2_zeroA),
             ('sig2_beta_T1', lambda x: x._params1.find_mean_sig2_beta()),
             ('sig2_beta_T2', lambda x: x._params2.find_mean_sig2_beta()),
             ('h2_T1', lambda x: x._params1.find_h2()),
             ('h2_T2', lambda x: x._params2.find_h2()),
             ('rho_zero', lambda x: x._rho_zero),
             ('rho_beta', lambda x: x._rho_beta),
             ('rg_sig2_factor', lambda x: x._rg_sig2_factor()),
             ('rg', lambda x: x._rg()),
             ('pi1', lambda x: x._params1.pi - x._pi12),
             ('pi2', lambda x: x._params2.pi - x._pi12),
             ('pi12', lambda x: x._pi12),
             ('pi1u', lambda x: x._params1.pi),
             ('pi2u', lambda x: x._params2.pi),
             ('dice', lambda x: (2 * x._pi12) / (x._params1.pi + x._params2.pi)),
             ('nc1', lambda x: x._params1._libbgmg.num_snp * (x._params1.pi - x._pi12)),
             ('nc2', lambda x: x._params1._libbgmg.num_snp * (x._params2.pi - x._pi12)),
             ('nc12', lambda x: x._params1._libbgmg.num_snp * x._pi12),
             ('nc1u', lambda x: x._params1._libbgmg.num_snp * x._params1.pi),
             ('nc2u', lambda x: x._params1._libbgmg.num_snp * x._params2.pi),
             ('nc1@p9', lambda x: NCKoef * x._params1._libbgmg.num_snp * (x._params1.pi - x._pi12)),
             ('nc2@p9', lambda x: NCKoef * x._params1._libbgmg.num_snp * (x._params2.pi - x._pi12)),
             ('nc12@p9', lambda x: NCKoef * x._params1._libbgmg.num_snp * x._pi12),
             ('totalpi', lambda x: (x._params1.pi+x._params2.pi-x._pi12)),
             ('totalnc', lambda x: x._params1._libbgmg.num_snp * (x._params1.pi+x._params2.pi-x._pi12)),
             ('totalnc@p9', lambda x: NCKoef * x._params1._libbgmg.num_snp * (x._params1.pi+x._params2.pi-x._pi12)),
             ('pi1_over_totalpi', lambda x: (x._params1.pi - x._pi12) / (x._params1.pi + x._params2.pi - x._pi12)),
             ('pi2_over_totalpi', lambda x: (x._params2.pi - x._pi12) / (x._params1.pi + x._params2.pi - x._pi12)),
             ('pi12_over_totalpi', lambda x: x._pi12 / (x._params1.pi + x._params2.pi - x._pi12)),
             ('pi12_over_pi1u', lambda x: x._pi12 / x._params1.pi),
             ('pi12_over_pi2u', lambda x: x._pi12 / x._params2.pi),
             ('pi1u_over_pi2u', lambda x: x._params1.pi / x._params2.pi),
             ('pi2u_over_pi1u', lambda x: x._params2.pi / x._params1.pi)]
              
    stats = [('mean', lambda x: np.mean(x)),
             ('median', lambda x: np.median(x)),
             ('std', lambda x: np.std(x)),
             ('lower', lambda x: np.percentile(x, 100.0 * (  alpha/2))),
             ('upper', lambda x: np.percentile(x, 100.0 * (1-alpha/2)))]

    return funcs, stats

def calc_qq_data(zvec, weights, hv_logp):
    assert(np.all(np.isfinite(zvec)))
    assert(np.all(np.isfinite(weights)))
    if np.sum(weights) == 0:
        retval = np.empty(hv_logp.shape)
        retval[:] = np.nan
        return retval

    # Step 0. calculate weights for all data poitns
    # Step 1. convert zvec to -log10(pvec)
    # Step 2. sort pval from large (1.0) to small (0.0), and take empirical cdf of this distribution
    # Step 3. interpolate (data_x, data_y) curve, as we don't need 10M data points on QQ plots
    data_weights = weights / np.sum(weights)                            # step 0
    data_y = -np.log10(2*scipy.stats.norm.cdf(-np.abs(zvec)))           # step 1
    si = np.argsort(data_y); data_y = data_y[si]                        # step 2
    data_x=-np.log10(np.flip(np.cumsum(np.flip(data_weights[si]))))     # step 2
    data_idx = np.not_equal(data_y, np.concatenate((data_y[1:], [np.inf])))
    data_logpvec = interp1d(data_y[data_idx], data_x[data_idx],         # step 3
                            bounds_error=False, fill_value=np.nan)(hv_logp)
    return data_logpvec

def calc_qq_model(zgrid, pdf, hv_z):
    model_cdf = np.cumsum(pdf) * (zgrid[1] - zgrid[0])
    model_cdf = 0.5 * (np.concatenate(([0.0], model_cdf[:-1])) + np.concatenate((model_cdf[:-1], [1.0])))
    model_logpvec = -np.log10(2*interp1d(-zgrid[zgrid<=0], model_cdf[zgrid<=0],
                                         bounds_error=False, fill_value=np.nan)(hv_z))
    return model_logpvec

def calc_qq_plot(libbgmg, params, trait_index, downsample_factor, mask=None, title=''):
    # mask can subset SNPs that are going into QQ curve, for example LDxMAF bin.
    if mask is None:
        mask = np.ones((libbgmg.num_tag, ), dtype=bool)

    zvec = libbgmg.get_zvec(trait_index)
    nvec = libbgmg.get_nvec(trait_index)
    weights = libbgmg.weights
    defvec = mask & (weights > 0) & np.isfinite(zvec) & np.isfinite(nvec)

    # Regular grid (vertical axis of the QQ plots)
    hv_z = np.linspace(0, 38, 1000)
    assert(np.all(np.isfinite(hv_z)))
    hv_logp = -np.log10(2*scipy.stats.norm.cdf(-hv_z))

    # Empirical (data) QQ plot
    data_logpvec = calc_qq_data(zvec[defvec], weights[defvec], hv_logp)

    # Estimated (model) QQ plots
    model_weights = downsample_weights(weights, downsample_factor, defvec, normalize=True)

    if np.sum(model_weights) == 0:
        model_logpvec = np.empty(hv_logp.shape)
        model_logpvec[:] = np.nan
    else:
        zgrid = np.arange(0, 38.0, 0.05, np.float32)
        try:
            libbgmg.weights = model_weights
            pdf = params.pdf(libbgmg, trait_index, zgrid)
        finally:
            libbgmg.weights = weights

        zgrid = np.concatenate((np.flip(-zgrid[1:]), zgrid))  # extend [0, 38] to [-38, 38]
        pdf = np.concatenate((np.flip(pdf[1:]), pdf))
        model_logpvec = calc_qq_model(zgrid, pdf, hv_z)

    return {'hv_logp': hv_logp,
            'data_logpvec': data_logpvec,
            'model_logpvec': model_logpvec,
            'n_snps': int(np.sum(defvec)),
            'sum_data_weights': float(np.sum(weights[defvec])),
            'title' : title}

def calc_bivariate_pdf(libbgmg, params, downsample_factor):
    weights = libbgmg.weights 
    defvec = weights > 0
    defvec = defvec & np.isfinite(libbgmg.zvec1) & np.isfinite(libbgmg.nvec1)
    defvec = defvec & np.isfinite(libbgmg.zvec2) & np.isfinite(libbgmg.nvec2)
    model_weights = downsample_weights(weights, downsample_factor, defvec, normalize=True)

    zgrid = np.arange(-25, 25.00001, 0.05)
    if np.sum(model_weights) == 0:
        raise(ValueError('np.sum(model_weights) == 0; consider reducing --downsample-factor ?'))

    try:
        libbgmg.weights = model_weights
        pdf = params.pdf(libbgmg, zgrid)
    finally:
        libbgmg.weights = weights

    return zgrid, pdf

def calc_bivariate_qq(libbgmg, zgrid, pdf):
    zgrid_fine = np.arange(-25, 25.00001, 0.005)   # project to a finer grid
    pdf_fine=scipy.interpolate.RectBivariateSpline(zgrid, zgrid, pdf)(zgrid_fine, zgrid_fine,grid=True)
    pthresh_vec = [1, 0.1, 0.01, 0.001]
    zthresh_vec = -scipy.stats.norm.ppf(np.array(pthresh_vec)/2)

    defvec = libbgmg.weights > 0
    defvec = defvec & np.isfinite(libbgmg.zvec1) & np.isfinite(libbgmg.nvec1)
    defvec = defvec & np.isfinite(libbgmg.zvec2) & np.isfinite(libbgmg.nvec2)

    zvec1=libbgmg.zvec1[defvec]
    zvec2=libbgmg.zvec2[defvec]
    weights = libbgmg.weights[defvec]

    # Regular grid (vertical axis of the QQ plots)
    hv_z = np.linspace(0, np.min([np.max(np.abs(np.concatenate((zvec1, zvec2)))), 38.0]), 1000)
    hv_logp = -np.log10(2*scipy.stats.norm.cdf(-hv_z))

    result = []
    for zthresh, pthresh in zip(zthresh_vec, pthresh_vec):
        mask = abs(zvec2)>=zthresh
        data_logpvec = calc_qq_data(zvec1[mask], weights[mask], hv_logp)

        pd_cond = np.sum(pdf_fine[abs(zgrid_fine) >= zthresh, :], axis=0)
        pd_cond = pd_cond / np.sum(pd_cond) / (zgrid_fine[1]-zgrid_fine[0])
        model_logpvec = calc_qq_model(zgrid_fine, pd_cond, hv_z)

        title = 'T1|T2|{}'.format(pthresh)
        result.append({'hv_logp': hv_logp,
                       'data_logpvec': data_logpvec,
                       'model_logpvec': model_logpvec,
                       'n_snps': int(np.sum(mask)),
                       'sum_data_weights': float(np.sum(weights[mask])),
                       'title' : title})

    for zthresh in zthresh_vec:
        mask = abs(zvec1)>=zthresh
        data_logpvec = calc_qq_data(zvec2[mask], weights[mask], hv_logp)

        pd_cond = np.sum(pdf_fine[:, abs(zgrid_fine) >= zthresh], axis=1)
        pd_cond = pd_cond / np.sum(pd_cond) / (zgrid_fine[1]-zgrid_fine[0])
        model_logpvec = calc_qq_model(zgrid_fine, pd_cond, hv_z)

        title = 'T2|T1|{}'.format(pthresh)
        result.append({'hv_logp': hv_logp,
                       'data_logpvec': data_logpvec,
                       'model_logpvec': model_logpvec,
                       'n_snps': int(np.sum(mask)),
                       'sum_data_weights': float(np.sum(weights[mask])),
                       'title' : title})
    return result

def downsample_weights(weights, downsample_factor, defvec, normalize=True):
    if np.sum(weights[defvec]) == 0:
        return weights.copy()
    idx = np.nonzero(defvec)[0]
    num_samples = np.random.binomial(n=len(idx), p=1.0/float(downsample_factor))
    idx_downsample = np.random.choice(idx, size=num_samples, replace=False)
    new_weights = np.zeros((len(weights), ), dtype=np.float32)
    new_weights[idx_downsample] = weights[idx_downsample]
    if np.sum(new_weights) == 0:
        new_weights = weights.copy()  # fall back to using the entire array
    if normalize:
        new_weights = new_weights / np.sum(new_weights)
    return new_weights

def calc_power_curve(libbgmg, params, trait_index, downsample_factor, nvec=None):
    power_nvec = np.power(10, np.arange(0, 9, 0.1)) if (nvec is None) else nvec

    zvec = libbgmg.get_zvec(trait_index)
    nvec = libbgmg.get_nvec(trait_index)
    weights = libbgmg.weights
    defvec = (weights > 0) & np.isfinite(zvec) & np.isfinite(nvec)

    model_weights = downsample_weights(weights, downsample_factor, defvec, normalize=True)
    if np.sum(model_weights) == 0:
        power_svec = np.empty(power_nvec.shape)
        power_svec[:] = np.nan
    else:
        try:
            libbgmg.weights = model_weights
            power_svec = params.power(libbgmg, trait_index, power_nvec)
        finally:
            libbgmg.weights = weights

    return {'nvec': power_nvec, 'svec': power_svec}

def _params_to_dict(params):
    if isinstance(params, BivariateParams):
        return {'pi12': params._pi12, 'rho_zero': params._rho_zero, 'rho_beta': params._rho_beta,
                'params1': _params_to_dict(params._params1), 'params2': _params_to_dict(params._params2),
                'type':'BivariateParams' }
    return gsa_mixer.utils._params_to_dict(params)

def _dict_to_params(p, libbgmg, args=None):
    if 'type' not in p:
        raise ValueError('Unrecognized type in dict_to_params()')
    if p['type'] == 'BivariateParams':
        return BivariateParams(params1=_dict_to_params(p['params1'], libbgmg, args),
                               params2=_dict_to_params(p['params2'], libbgmg, args),
                               pi12=p['pi12'], rho_beta=p['rho_beta'], rho_zero=p['rho_zero'])
    return gsa_mixer.utils._dict_to_params(p, libbgmg, args)
