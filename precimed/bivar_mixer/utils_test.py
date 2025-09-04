import numpy as np
import logging
import collections
from bivar_mixer.utils_obsolete import UnivariateParams_obsolete
from scipy.sparse import csr_matrix, coo_matrix, csc_matrix
import pytest

from .utils import UnivariateParametrization_natural_axis
from .utils import UnivariateParametrization_constPI
from .utils import UnivariateParametrization_constPI_constSIG2BETA
from gsa_mixer.utils import AnnotUnivariateParams
from .utils import BivariateParams
from .utils import _calculate_univariate_uncertainty
from .cli import apply_univariate_fit_sequence
from .cli import apply_bivariate_fit_sequence

def get_null_params(libbgmg):
    return AnnotUnivariateParams(libbgmg=libbgmg, s=0, l=0, sig2_zeroL=0)

def get_params(libbgmg):
    return get_null_params(libbgmg).copy(pi=0.02, sig2_beta=0.0005, sig2_zeroA=1.5)

def get_true_params(libbgmg):
    return get_null_params(libbgmg).copy(pi=0.1, sig2_beta=0.001, sig2_zeroA=1.2)

class LibBgmgMock(object):
    def __init__(self, num_snp, num_tag, make_weights=False):
        self.num_snp = num_snp
        self.num_tag = num_tag
        self._mafvec = np.random.uniform(0.01, 0.50, size=(num_snp,)).astype(np.float32)
        self._tldvec = np.random.uniform(1, 10, size=(num_snp,)).astype(np.float32)
        self.aux_gradient = np.random.normal(size=(self.num_snp + 2,)).astype(np.float32)
        self.weights = None
        if make_weights: self.weights = np.random.uniform(0, 1, size=(self.num_tag,)).astype(np.float32)
        self.cost = 0.0
        self.annomat = np.ones((num_snp, 1))
        self.annonames = ['base']
        self.genemat = np.ones((num_snp, 1))
        self.genenames = ['base']

    def set_option(self, option, value):
        pass

    def calc_unified_univariate_cost(self, trait, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL):
        if self.weights is None: return self.cost
        true_sig2_zeroA = get_true_params(self).sig2_zeroA
        true_pi = get_true_params(self).pi
        true_sig2beta = get_true_params(self).sig2_beta
        true_h2 = true_pi * true_sig2beta * self.num_snp
        pi_vec_flat = np.squeeze(np.asarray(pi_vec))
        sig2_beta_flat = np.squeeze(np.asarray(sig2_vec))
        h2_flat = np.multiply(pi_vec_flat, sig2_beta_flat)
        f1 = np.sum(np.power(pi_vec_flat - true_pi, 2))
        f2 = np.power(np.sum(h2_flat) - true_h2, 2)
        f4 = np.sum(np.power(sig2_beta_flat - true_sig2beta, 2))
        f3 = (sig2_zeroA - true_sig2_zeroA)**2
        cost = f1 + f2 + f3  + f4
        #print(pi_vec[0][0], sig2_vec[0][0], sig2_zeroA, cost)        
        return 10000*cost

    def calc_unified_univariate_cost_with_gradients(self, trait, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL):
        return self.cost, self.aux_gradient

    def calc_unified_bivariate_cost(self, pi_vec, sig2_beta, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL):
        return self.cost

    def log_message(self, message):
        logging.info('log_message({})'.format(message)) 

def check_optimize_result(optimize_result, expected_cost_df):
    assert(optimize_result['message'] in ['Optimization terminated successfully.', 'Not a bracketing interval.', 'Bracketing values (xa, xb, xc) do not fulfill this requirement: (f(xb) < f(xa)) and (f(xb) < f(xc))'])
    assert(optimize_result['nfev'] >= 4)
    assert(optimize_result['nit'] >= 1)
    assert(np.isfinite(optimize_result['AIC']))
    assert(np.isfinite(optimize_result['BIC']))
    if expected_cost_df is not None:
        assert(optimize_result['cost_df'] == expected_cost_df)

# py.test precimed/mixer/utils_test.py  -k test_univariate_fit_diffevo
def test_univariate_fit_diffevo():
    libbgmg_mock = LibBgmgMock(num_snp=10, num_tag=5, make_weights=True)
    args = collections.namedtuple('Args', ['diffevo_fast_repeats', 'seed', 'load_params_file'])._make((3, 123, ''))
    params, optimize_result_sequence = apply_univariate_fit_sequence(args, libbgmg_mock, ['diffevo-fast'], init_params=None, trait=1)
    print(params, optimize_result_sequence)
    assert(isinstance(params, AnnotUnivariateParams))
    check_optimize_result(optimize_result_sequence[0][1], expected_cost_df=3)

# py.test precimed/mixer/utils_test.py  -k test_univariate_fit_neldermead
def test_univariate_fit_neldermead():
    libbgmg_mock = LibBgmgMock(num_snp=10, num_tag=5, make_weights=True)
    args = collections.namedtuple('Args', ['diffevo_fast_repeats', 'seed', 'load_params_file'])._make((3, 123, ''))
    params, optimize_result_sequence = apply_univariate_fit_sequence(args, libbgmg_mock, ['neldermead-fast'], init_params=get_params(libbgmg_mock), trait=1)
    print(params, optimize_result_sequence)
    assert(isinstance(params, AnnotUnivariateParams))
    check_optimize_result(optimize_result_sequence[0][1], expected_cost_df=3)

# py.test precimed/mixer/utils_test.py  -k test_univariate_fit_inflation
def test_univariate_fit_inflation():
    libbgmg_mock = LibBgmgMock(num_snp=10, num_tag=5, make_weights=True)
    args = collections.namedtuple('Args', ['diffevo_fast_repeats', 'seed', 'load_params_file'])._make((3, 123, ''))
    params, optimize_result_sequence = apply_univariate_fit_sequence(args, libbgmg_mock, ['inflation'], init_params=get_params(libbgmg_mock), trait=1)
    print(params, optimize_result_sequence)
    assert(isinstance(params, AnnotUnivariateParams))
    check_optimize_result(optimize_result_sequence[0][1], expected_cost_df=1)

# py.test precimed/mixer/utils_test.py  -k test_univariate_fit_infinitesimal
def test_univariate_fit_infinitesimal():
    libbgmg_mock = LibBgmgMock(num_snp=10, num_tag=5, make_weights=True)
    args = collections.namedtuple('Args', ['diffevo_fast_repeats', 'seed', 'load_params_file'])._make((3, 123, ''))
    params, optimize_result_sequence = apply_univariate_fit_sequence(args, libbgmg_mock, ['infinitesimal'], init_params=get_params(libbgmg_mock), trait=1)
    print(params, optimize_result_sequence)
    assert(isinstance(params, AnnotUnivariateParams))
    check_optimize_result(optimize_result_sequence[0][1], expected_cost_df=2)

# py.test precimed/mixer/utils_test.py  -k test_univariate_uncertainty
def test_univariate_uncertainty():
    libbgmg_mock = LibBgmgMock(num_snp=10, num_tag=5, make_weights=True)    
    num_ci_samples = 5
    params = get_true_params(libbgmg_mock)
    NCKoef = params.find_nckoef()
    ci, ci_sample = _calculate_univariate_uncertainty(params, UnivariateParametrization_natural_axis(params, libbgmg_mock, trait=1), 0.95, num_samples=num_ci_samples, NCKoef=NCKoef)
    #print(ci, ci_sample)
    for k in ci.keys():
        for stat in ci[k].keys():
            assert(np.isfinite(ci[k][stat]))
    assert(len(ci_sample) == num_ci_samples)

# py.test precimed/mixer/utils_test.py  -k test_bivariate_fit_sequence
def test_bivariate_fit_sequence():
    libbgmg_mock = LibBgmgMock(num_snp=10, num_tag=5, make_weights=True)
    args = collections.namedtuple('Args', ['diffevo_fast_repeats', 'seed', 'analysis', 'load_params_file', 'trait1_params_file', 'trait2_params_file'])._make((3, 123, 'fit2', None, get_params(libbgmg_mock), get_params(libbgmg_mock)))
    params, _, _, optimize_result_sequence = apply_bivariate_fit_sequence(args, libbgmg_mock, ['diffevo-fast', 'neldermead-fast', 'brute1', 'brent1'])
    #print(params, optimize_result_sequence)
    assert(isinstance(params, BivariateParams))
    for optimize_result in optimize_result_sequence:
        check_optimize_result(optimize_result[1], expected_cost_df=None)
    assert(len(optimize_result_sequence) == 4)

    args = collections.namedtuple('Args', ['analysis', 'load_params_file'])._make(('test2', params))
    params, _, _, optimize_result_sequence = apply_bivariate_fit_sequence(args, libbgmg_mock, ['load', 'inflation'])
    print(params, optimize_result_sequence)
    for optimize_result in optimize_result_sequence:
        check_optimize_result(optimize_result[1], expected_cost_df=None)
    assert(len(optimize_result_sequence) == 3)  # inflation is fitted for sig2_zeroA(trait1), sig2_zeroA(trait2), and rho_zero

# py.test precimed/mixer/utils_test.py  -k test_parametrizations
@pytest.mark.parametrize("UnivariateParametrization,expected_vec",
    [(UnivariateParametrization_natural_axis, [0.20701416938432612, -9.210340371976182, -6.906754778648554]),
     (UnivariateParametrization_constPI, [-9.210340371976182, 0.20701416938432612]), 
     (UnivariateParametrization_constPI_constSIG2BETA, [0.20701416938432612])])
def test_parametrizations(UnivariateParametrization, expected_vec):
    params = UnivariateParams_obsolete(0.001, 1e-4, 1.23)
    pp=UnivariateParametrization(params, None, 1)
    vec=pp.params_to_vec(params)
    assert(np.all(np.isclose(vec, expected_vec)))
    params2 = pp.vec_to_params(vec)
    assert(np.all(np.isclose([params.pi, params.sig2_beta, params.sig2_zeroA],
                             [params2.pi, params2.sig2_beta, params2.sig2_zeroA])))
