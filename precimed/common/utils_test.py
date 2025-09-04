import numpy as np

from .utils import _log_exp_converter
from .utils import _logit_logistic_converter
from .utils import _arctanh_tanh_converter

# py.test precimed/common/utils_test.py  -k test_converters
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
