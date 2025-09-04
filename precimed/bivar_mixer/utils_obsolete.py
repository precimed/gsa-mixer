import numpy as np

from common.utils import _arctanh_tanh_converter, _log_exp_converter

class UnivariateParams_obsolete(object):
    def __init__(self, pi, sig2_beta, sig2_zeroA):
        self._pi = pi
        self._sig2_beta = sig2_beta
        self._sig2_zeroA = sig2_zeroA
        self._validate()

    def copy(self):
        return UnivariateParams_obsolete(self._pi, self._sig2_beta, self._sig2_zeroA)

    @property
    def pi(self):
        return self._pi

    @property
    def sig2_beta(self):
        return self._sig2_beta

    @property
    def sig2_zeroA(self):
        return self._sig2_zeroA

    def _validate(self):
        for val in [self._pi, self._sig2_beta, self._sig2_zeroA]:
            assert np.isscalar(val)
            assert np.greater_equal(val, 0)
        assert np.less_equal(self._pi, 1)
    
    def __str__(self):
        description = []
        for attr_name in '_pi', '_sig2_beta', '_sig2_zeroA':
            try:
                attr_value = getattr(self, attr_name)
                description.append('{}: {}'.format(attr_name, attr_value))
            except RuntimeError:
                pass
        return 'UnivariateParams({})'.format(', '.join(description))
    __repr__ = __str__

    def find_pi_mat(self, num_snp):
        return self._pi * np.ones(shape=(num_snp, 1), dtype=np.float32)

    def find_sig2_mat(self, num_snp):
        return self._sig2_beta * np.ones(shape=(num_snp, 1), dtype=np.float32)

    def cost(self, lib, trait):
        value = lib.calc_unified_univariate_cost(trait, self.find_pi_mat(lib.num_snp), self.find_sig2_mat(lib.num_snp), 
                                                 sig2_zeroA=self._sig2_zeroA, sig2_zeroC=1, sig2_zeroL=0)
        return value if np.isfinite(value) else 1e100

    def aux(self, lib, trait):
        return lib.calc_unified_univariate_aux(trait, self.find_pi_mat(lib.num_snp), self.find_sig2_mat(lib.num_snp), 
                                                 sig2_zeroA=self._sig2_zeroA, sig2_zeroC=1, sig2_zeroL=0)

    def pdf(self, lib, trait, zgrid):
        return lib.calc_unified_univariate_pdf(trait, self.find_pi_mat(lib.num_snp), self.find_sig2_mat(lib.num_snp),
                                               sig2_zeroA=self._sig2_zeroA, sig2_zeroC=1, sig2_zeroL=0, zgrid=zgrid)

    def power(self, lib, trait, ngrid, zthresh=5.45):
        svec_num, svec_denom = lib.calc_unified_univariate_power(trait, self.find_pi_mat(lib.num_snp), self.find_sig2_mat(lib.num_snp),
                                                 sig2_zeroA=self._sig2_zeroA, sig2_zeroC=1, sig2_zeroL=0, zthresh=zthresh, ngrid=ngrid)
        return np.divide(svec_num, svec_denom)

# Unconstrained parametrization with "independent axis", i.e.
#   x1 = log(sig2_zeroA)
#   x2 = log(atanh(pi)) + log(sig2_beta)
#   x3 = log(atanh(pi)) - log(sig2_beta)
# The reason for atanh(pi) in the formulas is to make sure that the inverse transform always give a valid pi (e.g. between zero to one)
# atanh is particularly helpful because atanh(x)~x for small x, and we expect pi to be small.
class UnivariateParametrization(object):
    def __init__(self, params, lib, trait):
        self._params = params.copy()
        self._lib = lib
        self._trait = trait

    def params_to_vec(self, params):
        arctanh_pi = _arctanh_tanh_converter(params.pi, invflag=False)
        return [_log_exp_converter(params._sig2_zeroA, invflag=False),
                _log_exp_converter(arctanh_pi, invflag=False) + _log_exp_converter(params._sig2_beta, invflag=False),
                _log_exp_converter(arctanh_pi, invflag=False) - _log_exp_converter(params._sig2_beta, invflag=False)]
        
    def vec_to_params(self, vec):
        self._params._sig2_zeroA=_log_exp_converter(vec[0], invflag=True)
        self._params._sig2_beta=_log_exp_converter((vec[1] - vec[2]) / 2.0, invflag=True)
        self._params._pi= _arctanh_tanh_converter(_log_exp_converter((vec[1] + vec[2]) / 2.0, invflag=True), invflag=True)
        return self._params
    
    def calc_cost(self, vec):
        return self.vec_to_params(vec).cost(self._lib, self._trait)

