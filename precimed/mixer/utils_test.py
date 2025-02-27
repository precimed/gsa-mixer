import numpy as np
import logging
import collections
from precimed.mixer.utils_obsolete import UnivariateParams_obsolete
from scipy.sparse import csr_matrix, coo_matrix, csc_matrix
import pytest

from .utils import UnivariateParametrization_natural_axis
from .utils import UnivariateParametrization_constPI
from .utils import UnivariateParametrization_constPI_constSIG2BETA
from .utils import AnnotUnivariateParams
from .utils import AnnotUnivariateParametrization
from .utils import BivariateParams
from .utils import _calculate_univariate_uncertainty
from .cli import apply_univariate_fit_sequence
from .cli import apply_bivariate_fit_sequence

from .utils import _log_exp_converter
from .utils import _logit_logistic_converter
from .utils import _arctanh_tanh_converter

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

def make_test_objects(num_annot, num_gene, num_snp, num_tag, seed=123):
    np.random.seed(seed)
    libbgmg = LibBgmgMock(num_snp, num_tag)
    libbgmg.annomat = np.random.binomial(1, 0.5, (num_snp, num_annot))
    libbgmg.annonames = ['annot{}'.format(i + 1) for i in range(num_annot)]
    libbgmg.genemat = np.random.binomial(1, 0.5, (num_snp, num_gene)) if (num_gene > 1) else np.ones((num_snp, 1))
    libbgmg.genenames = ['gene{}'.format(i + 1) for i in range(num_gene)] if (num_gene > 1) else ['base']
    return libbgmg

def assert_params(params, params2):
    assert(np.isclose(params2._pi, params._pi))
    assert(np.isclose(params2._sig2_beta, params._sig2_beta))
    assert(np.all(np.isclose(params2._sig2_annot, params._sig2_annot)))
    assert(np.isclose(params2._s, params._s))
    assert(np.isclose(params2._l, params._l))
    assert(np.isclose(params2._sig2_zeroA, params._sig2_zeroA))
    assert(np.isclose(params2._sig2_zeroL, params._sig2_zeroL))

# py.test precimed/mixer/utils_test.py  -k test_AnnotUnivariateParams
def test_AnnotUnivariateParams():
    num_annot = 4; num_gene = 1; num_snp = 20; num_tag = 20; trait_index = 1
    libbgmg = make_test_objects(num_annot, num_gene, num_snp, num_tag)
    libbgmg.annomat = csr_matrix(libbgmg.annomat.astype(np.float32))
    libbgmg.genemat = csr_matrix(libbgmg.genemat.astype(np.float32))
    params = AnnotUnivariateParams(pi=1e-3, sig2_beta=1, sig2_annot=np.random.uniform(0.01, 0.1, size=(len(libbgmg.annonames),)), sig2_gene=[1], s=-0.5, l=-0.25, sig2_zeroA=1.23, sig2_zeroL=1e-5, libbgmg=libbgmg)
    constr = AnnotUnivariateParams(pi=1e-3, sig2_beta=1, sig2_annot=[None for x in libbgmg.annonames], sig2_gene=[1], s=None, l=None, sig2_zeroA=None, sig2_zeroL=None, libbgmg=libbgmg)
    parametrization = AnnotUnivariateParametrization(trait_index, constr)

    sig2_mat = params.find_sig2_mat()
    assert(sig2_mat.shape[0] == num_snp)
    assert(sig2_mat.shape[1] == 1)

    cost = params.cost(libbgmg, trait_index)
    assert(cost==0.0)
    
    cost, aux = params.cost_with_gradients(libbgmg, trait_index)
    assert(len(aux) == (num_snp + 2))
    assert(cost==0.0)

    vec = parametrization.params_to_vec(params)
    assert((len(vec) == (num_annot + 4)) and np.all(np.isfinite(vec)))

    params2 = parametrization.vec_to_params(vec)
    assert_params(params, params2)

    cost, gradients = parametrization.calc_cost_with_gradients(vec, [0.1 for x in vec], verbose=True, force_numeric=False)
    expected_gradients = [2.79372632e+00, -1.08375253e+02,  2.57886913e+02, -7.73860615e+00,  2.79745433e+02, -1.64321041e+02, 1.83964229e+00,  1.06939267e-05]
    assert(np.all(np.isclose(gradients.flatten(), np.array(expected_gradients).astype(np.float32))))

    expected_h2 = [[0.3709366,  0.47575423, 0.3714563,  0.41520426]] 
    expected_snps = [[11, 11,  8, 12]]
    expected_enrich =  [[1,        1.2825756, 1.0014011, 1.1193402]]

    annot_h2 = params.find_annot_h2()
    annot_snps = params.find_annot_snps()
    annot_enrich = params.find_annot_enrich()
    
    assert(np.all(np.isclose(annot_h2, np.array(expected_h2).astype(np.float32))))
    assert(np.all(np.isclose(annot_snps, np.array(expected_snps).astype(np.float32))))
    assert(np.all(np.isclose(annot_enrich, np.array(expected_enrich).astype(np.float32))))

# py.test precimed/mixer/utils_test.py  -k test_AnnotUnivariateParamsWithGenemat
def test_AnnotUnivariateParamsWithGenemat():
    num_annot = 4; num_gene = 6; num_snp = 20; num_tag = 20; trait_index = 1
    libbgmg = make_test_objects(num_annot, num_gene, num_snp, num_tag)
    libbgmg.annomat = csr_matrix(libbgmg.annomat.astype(np.float32))
    libbgmg.genemat = csr_matrix(libbgmg.genemat.astype(np.float32))
    params = AnnotUnivariateParams(pi=1e-3, sig2_beta=1, sig2_annot=np.random.uniform(0.01, 0.1, size=(len(libbgmg.annonames),)), sig2_gene=np.random.uniform(0.01, 0.1, size=(len(libbgmg.genenames),)), s=-0.5, l=-0.25, sig2_zeroA=1.23, sig2_zeroL=1e-5, libbgmg=libbgmg)
    constr = AnnotUnivariateParams(pi=1e-3, sig2_beta=1, sig2_annot=[None for x in libbgmg.annonames], sig2_gene=[None for x in libbgmg.genenames], s=None, l=None, sig2_zeroA=None, sig2_zeroL=None, libbgmg=libbgmg)
    parametrization = AnnotUnivariateParametrization(trait_index, constr)

    sig2_mat = params.find_sig2_mat()
    assert(sig2_mat.shape[0] == num_snp)
    assert(sig2_mat.shape[1] == 1)

    cost = params.cost(libbgmg, trait_index)
    assert(cost==0.0)
    
    cost, aux = params.cost_with_gradients(libbgmg, trait_index)
    assert(len(aux) == (num_snp + 2))
    assert(cost==0.0)

    vec = parametrization.params_to_vec(params)
    assert((len(vec) == (num_annot + num_gene + 4)) and np.all(np.isfinite(vec)))

    params2 = parametrization.vec_to_params(vec)
    assert_params(params, params2)

    cost, gradients = parametrization.calc_cost_with_gradients(vec, [0.1 for x in vec], verbose=True, force_numeric=False)
    expected_gradients = [-1.4731736e+01, -8.2612171e+00,  5.5556526e+00, -1.1741324e+01,  3.3541510e+00, -1.0850573e+01, -1.1283661e+01,  6.4700727e+00,   1.2052662e+00, -1.8073883e+01,  3.3172977e+01, -1.9671848e+01, 1.8396423e+00,  1.0693927e-05]
    assert(np.all(np.isclose(gradients.flatten(), np.array(expected_gradients).astype(np.float32))))

    expected_h2 = [[0.02835298, 0.03408064, 0.02363124, 0.03536456]] 
    expected_snps = [[11, 11,  8, 12]]
    expected_enrich =  [[1. ,       1.2020128, 0.8334659, 1.2472961]]
    annot_h2 = params.find_annot_h2(params._annomat)
    annot_snps = params.find_annot_snps(params._annomat)
    annot_enrich = params.find_annot_enrich(params._annomat)
    #print(annot_h2);    print(annot_snps);    print(annot_enrich)
    assert(np.all(np.isclose(annot_h2, np.array(expected_h2).astype(np.float32))))
    assert(np.all(np.isclose(annot_snps, np.array(expected_snps).astype(np.float32))))
    assert(np.all(np.isclose(annot_enrich, np.array(expected_enrich).astype(np.float32))))

    expected_h2 = [[0.03100209, 0.02786526, 0.03020591, 0.02881684, 0.02202258, 0.0371193]] 
    expected_snps = [[10, 9, 11,  9, 8, 13]]
    expected_enrich =  [[1,       0.8988186,  0.9743185,  0.92951286, 0.7103579,  1.197316]]
    gene_h2 = params.find_annot_h2(params._genemat)
    gene_snps = params.find_annot_snps(params._genemat)
    gene_enrich = params.find_annot_enrich(params._genemat)
    #print(gene_h2);    print(gene_snps);    print(gene_enrich)
    assert(np.all(np.isclose(gene_h2, np.array(expected_h2).astype(np.float32))))
    assert(np.all(np.isclose(gene_snps, np.array(expected_snps).astype(np.float32))))
    assert(np.all(np.isclose(gene_enrich, np.array(expected_enrich).astype(np.float32))))

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

# py.test precimed/mixer/utils_test.py  -k test_converters
def test_converters():
    epsval = np.finfo(float).eps
    minval = np.finfo(float).min
    maxval = np.finfo(float).max

    def isclose(a, b, rel_tol=1e-09, abs_tol=2*epsval):
        return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

    tests = [(_log_exp_converter,        0.0, maxval, [0.00001, 0.001, 1.2, 1000, 10000]), 
             (_logit_logistic_converter, 0.0, 1.0,   [0.00001, 0.001, 0.1, 0.9, 0.999, 0.99999]),
             (_arctanh_tanh_converter,    -1.0, 1.0,  [-0.9999, -0.991, -0.9, 0.9, 0.999, 0.99999])]

    for func, limLow, limHigh, values in tests:
        for index, val in enumerate([limLow] + values + [limHigh]):
            #print('{}: {} -> {} -> {}'.format(func, val, func(val, False), func(func(val, False), True)))
            assert np.isfinite(func(val, False)), 'Error in {}({})'.format(func, val)
            assert isclose(val, func(func(val, False), True)), 'Error in {}({})'.format(func, val)
            isLimit = (index==0) or (index==(len(values)+ 1))
            assert abs(func(val, False)) < (1000 if isLimit else 20), '{}({}) is too large ({})'.format(func, val, func(val, False))
            assert abs(func(val, False)) > (0.001 if isLimit else 0.1), '{}({}) is too small ({})'.format(func, val, func(val, False))
        assert isclose(func(-10000, True), limLow), '{} vs {}'.format(func(-10000, True), limLow)
        assert isclose(func( 10000, True), limHigh), '{} vs {}'.format(func(10000, True), limHigh)

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
