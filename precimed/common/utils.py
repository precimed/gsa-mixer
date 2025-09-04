import json
import numpy as np

epsval = np.finfo(float).eps
minval = np.finfo(float).min
maxval = np.finfo(float).max

def _log_bounded(x): # transforms, roughtly speaking, [0.00001, 100000] to [-10, 10]
    if not np.isfinite(x): return x
    if x<epsval: x=epsval
    if x>maxval: x=maxval
    y = np.log(x)
    return y

def _exp_bounded(x): # transforms, roughtly speaking, [-10, 10] to [0.0001, 100000]
    if not np.isfinite(x): return x
    y = np.exp(x)
    if y<epsval: y=epsval
    if y>maxval: y=maxval
    return y

def _logit_bounded(x): # transforms, roughly speaking, [0.00001, 0.99999] to [-10, 10]
    if not np.isfinite(x): return x
    if x<epsval: x=epsval
    if x>(1-epsval): x=1-epsval
    y = _log_bounded(x / (1-x))
    return y
    
def _logistic_bounded(x): # transforms, roughly speaking, [-10, 10] to [0.00001, 0.99999]
    if not np.isfinite(x): return x
    y = _exp_bounded(x) / (1 + _exp_bounded(x))
    if y<epsval: y=epsval
    if y>(1-epsval): y=1-epsval
    return y
    
# converter that maps [0, +inf] domain to [-inf, inf], and back if infvlag=True
def _log_exp_converter(x, invflag=False):
    return _exp_bounded(x) if invflag else _log_bounded(x)

# converter that maps [0, 1] domain to [-inf, inf], and back if infvlag=True
def _logit_logistic_converter(x, invflag=False):
    return _logistic_bounded(x) if invflag else _logit_bounded(x)

# converter that maps [-1, 1] domain to [-inf, inf], and back if infvlag=True
# this has an additional property that arctanh(X) ~ X for small values of X.
def _arctanh_tanh_converter(x, invflag=False):
    return (2.0 * _logistic_bounded(2.0 * x) - 1.0) if invflag else 0.5*_logit_bounded(0.5*x + 0.5)

def make_one_zeros_vec(length):
        retval = np.zeros((length, 1))
        retval[0] = 1
        return retval

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if callable(obj):
            return str(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.float32):
            return np.float64(obj)
        if isinstance(obj, np.uint32):
            return str(obj)
        if isinstance(obj, np.int64):
            return str(obj)
        return json.JSONEncoder.default(self, obj)
