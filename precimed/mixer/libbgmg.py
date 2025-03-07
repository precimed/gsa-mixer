import os
import sys
import ctypes
import six
import logging
import numpy as np

def _p2n(arg_value):  # python2native
    if arg_value==None:
        return ctypes.create_string_buffer("".encode('utf-8'))
    if isinstance(arg_value, str):
        return ctypes.create_string_buffer(arg_value.encode('utf-8'))
    return arg_value

def _n2p(res_value):  # native2python
    if isinstance(res_value, bytes):
        if six.PY3:
            return res_value.decode('utf-8')
    return res_value

class LibBgmg(object):
    def __init__(self, lib_name=None, context_id=0, init_log=None, dispose=False):
        self._context_id = context_id
        self.cdll, self._lib_name = self._load_cdll(lib_name)
        logging.info('__init__(lib_name={}, context_id={})'.format(self._lib_name, context_id))

        # dark magic - https://github.com/numpy/numpy/issues/6239
        _float32_pointer_type_base = np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS')
        _int32_pointer_type_base = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS')
        _int64_pointer_type_base = np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags='C_CONTIGUOUS')
        def _from_param_float32(cls, obj):
            return obj if (obj is None) else _float32_pointer_type_base.from_param(obj)
        def _from_param_int32(cls, obj):
            return obj if (obj is None) else _int32_pointer_type_base.from_param(obj)
        def _from_param_int64(cls, obj):
            return obj if (obj is None) else _int64_pointer_type_base.from_param(obj)
        float32_pointer_type = type('float32_pointer_type', (_float32_pointer_type_base,), {'from_param': classmethod(_from_param_float32)} )
        int32_pointer_type = type('int32_pointer_type', (_int32_pointer_type_base,), {'from_param': classmethod(_from_param_int32)} )
        int64_pointer_type = type('int64_pointer_type', (_int64_pointer_type_base,), {'from_param': classmethod(_from_param_int64)} )

        # set function signatures ('restype' and 'argtype') for all functions that involve non-integer types
        # (pointers, floats, doubles, etc - either as input or as output)
        self.cdll.bgmg_get_last_error.restype = ctypes.c_char_p
        self.cdll.bgmg_status.restype = ctypes.c_char_p
        self.cdll.bgmg_set_tag_indices.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, int32_pointer_type]
        self.cdll.bgmg_retrieve_tag_indices.argtypes = [ctypes.c_int, ctypes.c_int, int32_pointer_type]
        self.cdll.bgmg_set_mafvec.argtypes = [ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_retrieve_mafvec.argtypes = [ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_set_weights.argtypes = [ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_retrieve_weights.argtypes = [ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_set_chrnumvec.argtypes = [ctypes.c_int, ctypes.c_int, int32_pointer_type]
        self.cdll.bgmg_retrieve_chrnumvec.argtypes = [ctypes.c_int, ctypes.c_int, int32_pointer_type]
        self.cdll.bgmg_set_zvec.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_retrieve_zvec.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_set_nvec.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_retrieve_nvec.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_set_causalbetavec.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_retrieve_causalbetavec.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_retrieve_fixed_effect_delta.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_set_option.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_double]
        self.cdll.bgmg_set_ld_r2_coo_from_file.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_char_p]
        self.cdll.bgmg_set_ld_r2_csr.argtypes = [ctypes.c_int, ctypes.c_int]
        self.cdll.bgmg_set_weights_hardprune.argtypes = [ctypes.c_int, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_int]
        self.cdll.bgmg_set_weights_randprune.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_float, ctypes.c_float, ctypes.c_int, ctypes.c_char_p, ctypes.c_char_p]
        self.cdll.bgmg_perform_ld_clump.argtypes = [ctypes.c_int, ctypes.c_float, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_retrieve_ld_sum_r2.argtypes = [ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_retrieve_ld_sum_r4.argtypes = [ctypes.c_int, ctypes.c_int, float32_pointer_type]
        self.cdll.bgmg_num_ld_r_tag.argtypes = [ctypes.c_int, ctypes.c_int]
        self.cdll.bgmg_retrieve_ld_r_tag.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, int32_pointer_type, float32_pointer_type]
        self.cdll.bgmg_num_ld_r_chr.argtypes = [ctypes.c_int, ctypes.c_int]
        self.cdll.bgmg_retrieve_ld_r_chr.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_longlong, int32_pointer_type, int32_pointer_type, float32_pointer_type]
        self.cdll.bgmg_num_ld_r_tag_range.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int]
        self.cdll.bgmg_retrieve_ld_r_tag_range.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_longlong, int32_pointer_type, int32_pointer_type, float32_pointer_type]
        self.cdll.bgmg_calc_unified_univariate_cost.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, float32_pointer_type, float32_pointer_type, ctypes.c_float, ctypes.c_float, ctypes.c_float, float32_pointer_type]
        self.cdll.bgmg_calc_unified_univariate_cost.restype = ctypes.c_double
        self.cdll.bgmg_calc_unified_univariate_pdf.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, float32_pointer_type, float32_pointer_type, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_int, float32_pointer_type, float32_pointer_type]
        self.cdll.bgmg_calc_unified_univariate_power.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, float32_pointer_type, float32_pointer_type, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_int, float32_pointer_type, float32_pointer_type, float32_pointer_type]
        self.cdll.bgmg_calc_unified_univariate_delta_posterior.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, float32_pointer_type, float32_pointer_type, ctypes.c_float, ctypes.c_float, ctypes.c_float, ctypes.c_int, float32_pointer_type, float32_pointer_type, float32_pointer_type]

        self.cdll.bgmg_calc_unified_bivariate_cost.argtypes = [ctypes.c_int,            #int context_id
                                                               ctypes.c_int,            #int num_snp
                                                               float32_pointer_type,    #float* pi_vec
                                                               float32_pointer_type,    #float* sig2_vec
                                                               float32_pointer_type,    #float* rho_vec
                                                               float32_pointer_type,    #float* sig2_zeroA
                                                               float32_pointer_type,    #float* sig2_zeroC
                                                               float32_pointer_type,    #float* sig2_zeroL
                                                               ctypes.c_float,          #float rho_zeroA
                                                               ctypes.c_float,          #float rho_zeroL
                                                               float32_pointer_type]    #float* aux
        self.cdll.bgmg_calc_unified_bivariate_cost.restype = ctypes.c_double
        self.cdll.bgmg_calc_unified_bivariate_pdf.argtypes = [ctypes.c_int, ctypes.c_int, float32_pointer_type, float32_pointer_type, float32_pointer_type, float32_pointer_type, float32_pointer_type, float32_pointer_type, ctypes.c_float, ctypes.c_float, ctypes.c_int, float32_pointer_type, float32_pointer_type, float32_pointer_type]
        self.cdll.bgmg_calc_unified_bivariate_delta_posterior.argtypes = [ctypes.c_int, ctypes.c_int, float32_pointer_type, float32_pointer_type, float32_pointer_type, float32_pointer_type, float32_pointer_type, float32_pointer_type, ctypes.c_float, ctypes.c_float, ctypes.c_int, float32_pointer_type, float32_pointer_type, float32_pointer_type, float32_pointer_type, float32_pointer_type, float32_pointer_type]

        self.cdll.bgmg_calc_ld_matrix.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_float]

        self.cdll.bgmg_request_annotate_intervals.argtypes = [ctypes.c_int, int64_pointer_type, ctypes.c_int, int64_pointer_type, int64_pointer_type, int32_pointer_type]
        self.cdll.bgmg_request_annotate_intervals.restype = ctypes.c_int64
        self.cdll.bgmg_copy_requested_object.argtypes = [ctypes.c_int64, ctypes.c_char_p]

        if init_log: self.init_log(init_log)
        if dispose: self.dispose()

    def calc_ld_matrix(self, bfile, outfile, r2min, ldscore_r2min, ld_window, ld_window_kb):
        self.cdll.bgmg_calc_ld_matrix(_p2n(bfile), _p2n(outfile), r2min, ldscore_r2min, ld_window, np.float32(ld_window_kb))

    def get_last_error(self):
        return _n2p(self.cdll.bgmg_get_last_error())

    def init_log(self, file):
        logging.info('init_log({})'.format(file)) 
        self.cdll.bgmg_init_log(_p2n(file))

    def log_message(self, message):
        logging.info('log_message({})'.format(message)) 
        self.cdll.bgmg_log_message(_p2n(message))

    def dispose(self):
        return self._check_error(self.cdll.bgmg_dispose(self._context_id))

    def save(self, file):
        self.cdll.bgmg_save_context(self._context_id, _p2n(file))

    def load(self, file):
        self.cdll.bgmg_load_context(self._context_id, _p2n(file))

    @property
    def status(self):
        return _n2p(self.cdll.bgmg_status(self._context_id))

    def init(self, bim_file, frq_file, chr_labels, trait1_file, trait2_file, exclude, extract, exclude_ranges):
        chr_labels_val = chr_labels if isinstance(chr_labels, str) else ' '.join([str(x) for x in chr_labels])
        return self._check_error(self.cdll.bgmg_init(
            self._context_id, _p2n(bim_file), _p2n(frq_file), _p2n(chr_labels_val), _p2n(trait1_file), _p2n(trait2_file), _p2n(exclude), _p2n(extract), _p2n(exclude_ranges)))

    def read_trait_file(self, trait_index, trait_file, exclude, extract, exclude_ranges):
        return self._check_error(self.cdll.bgmg_read_trait_file(self._context_id, trait_index, _p2n(trait_file), _p2n(exclude), _p2n(extract), _p2n(exclude_ranges)))

    def set_option(self, option, value):
        if value is None: return None
        return self._check_error(self.cdll.bgmg_set_option(self._context_id, _p2n(option), value))

    def convert_plink_ld(self, plink_ld_gz, plink_ld_bin):
        return self._check_error(self.cdll.bgmg_convert_plink_ld(self._context_id, _p2n(plink_ld_gz), _p2n(plink_ld_bin)))

    def set_ld_r2_coo_from_file(self, chr_label, filename):
        return self._check_error(self.cdll.bgmg_set_ld_r2_coo_from_file(self._context_id, chr_label, _p2n(filename)))

    def set_ld_r2_csr(self, chr_label=-1):  # -1 means to finalize all chromosomes
        return self._check_error(self.cdll.bgmg_set_ld_r2_csr(self._context_id, chr_label))

    def set_weights_hardprune(self, subset, r2, maf, use_w_ld):
        return self._check_error(self.cdll.bgmg_set_weights_hardprune(self._context_id, subset, r2, maf, int(use_w_ld)))

    def set_weights_randprune(self, n, r2, maf, use_w_ld, exclude="", extract=""):
        return self._check_error(self.cdll.bgmg_set_weights_randprune(self._context_id, n, r2, maf, int(use_w_ld), _p2n(exclude), _p2n(extract)))

    def perform_ld_clump(self, r2, buffer): # buffer must have a length equal to 'num_tag'
        buffer_data = (buffer if isinstance(buffer, np.ndarray) else np.array(buffer)).astype(np.float32)
        self._check_error(self.cdll.bgmg_perform_ld_clump(self._context_id, r2, np.size(buffer), buffer_data))
        return buffer_data

    @property
    def num_tag(self):
        return self._check_error(self.cdll.bgmg_get_num_tag(self._context_id))
   
    @property
    def num_snp(self):
        return self._check_error(self.cdll.bgmg_get_num_snp(self._context_id))

    @property
    def defvec(self):
        numpy_ndarray = np.zeros(shape=(self.num_tag,), dtype=np.int32)
        self._check_error(self.cdll.bgmg_retrieve_tag_indices(self._context_id, np.size(numpy_ndarray), numpy_ndarray))
        mask_array = np.zeros(shape=(self.num_snp,), dtype=bool)
        mask_array[numpy_ndarray] = 1
        return mask_array

    @defvec.setter
    def defvec(self, val):
        indices = np.flatnonzero(np.array(val, dtype=bool)).astype(np.int32)
        self._check_error(self.cdll.bgmg_set_tag_indices(self._context_id, np.size(val), np.size(indices), indices))

    @property
    def mafvec(self):
        return self._get_vec_impl(self.cdll.bgmg_retrieve_mafvec, np.float32, self.num_snp, trait=None)

    @mafvec.setter
    def mafvec(self, val):
        self._set_vec_impl(self.cdll.bgmg_set_mafvec, np.float32, val, trait=None)

    @property
    def weights(self):
        return self._get_vec_impl(self.cdll.bgmg_retrieve_weights, np.float32, self.num_tag, trait=None)

    @weights.setter
    def weights(self, val):
        self._set_vec_impl(self.cdll.bgmg_set_weights, np.float32, val, trait=None)

    @property
    def chrnumvec(self):
        return self._get_vec_impl(self.cdll.bgmg_retrieve_chrnumvec, np.int32, self.num_snp, trait=None)

    @chrnumvec.setter
    def chrnumvec(self, val):  
        self._set_vec_impl(self.cdll.bgmg_set_chrnumvec, np.int32, val, trait=None)

    def set_zvec(self, val, trait):
        self._set_vec_impl(self.cdll.bgmg_set_zvec, np.float32, val, trait=trait)

    def get_zvec(self, trait):
        return self._get_vec_impl(self.cdll.bgmg_retrieve_zvec, np.float32, self.num_tag, trait=trait)

    def set_nvec(self, val, trait):
        self._set_vec_impl(self.cdll.bgmg_set_nvec, np.float32, val, trait=trait)

    def get_nvec(self, trait):
        return self._get_vec_impl(self.cdll.bgmg_retrieve_nvec, np.float32, self.num_tag, trait=trait)

    def set_causalbetavec(self, val, trait):
        self._set_vec_impl(self.cdll.bgmg_set_causalbetavec, np.float32, val, trait=trait)

    def get_causalbetavec(self, trait):
        return self._get_vec_impl(self.cdll.bgmg_retrieve_causalbetavec, np.float32, self.num_snp, trait=trait)

    def get_fixedeffectdelta(self, trait):
        return self._get_vec_impl(self.cdll.bgmg_retrieve_fixed_effect_delta, np.float32, self.num_snp, trait=trait)

    @property
    def zvec1(self):
        return self.get_zvec(trait=1)

    @zvec1.setter
    def zvec1(self, val):
        self.set_zvec(val, trait=1)

    @property
    def zvec2(self):
        return self.get_zvec(trait=2)

    @zvec2.setter
    def zvec2(self, val):
        self.set_zvec(val, trait=2)

    @property
    def nvec1(self):
        return self.get_nvec(trait=1)

    @nvec1.setter
    def nvec1(self, val):
        self.set_nvec(val, trait=1)

    @property
    def nvec2(self):
        return self.get_nvec(trait=2)

    @nvec2.setter
    def nvec2(self, val):
        self.set_nvec(val, trait=2)

    @property
    def causalbetavec1(self):
        return self.get_causalbetavec(trait=1)

    @causalbetavec1.setter
    def causalbetavec1(self, val):
        self.set_causalbetavec(val, trait=1)

    @property
    def causalbetavec2(self):
        return self.get_causalbetavec(trait=2)

    @causalbetavec2.setter
    def causalbetavec2(self, val):
        self.set_causalbetavec(val, trait=2)

    @property
    def ld_sum_r2(self):
        self.set_option("retrieve_ld_sum_type", 0)
        above_r2min = self._get_vec_impl(self.cdll.bgmg_retrieve_ld_sum_r2, np.float32, self.num_snp, trait=None)
        self.set_option("retrieve_ld_sum_type", 1)
        below_r2min = self._get_vec_impl(self.cdll.bgmg_retrieve_ld_sum_r2, np.float32, self.num_snp, trait=None)
        return above_r2min + below_r2min

    @property
    def ld_sum_r4(self):
        return self._get_vec_impl(self.cdll.bgmg_retrieve_ld_sum_r4, np.float32, self.num_snp, trait=None)

    # return (snp, r) tuple, representing LD structure of a given tag SNP
    # snp and r are vectors of the same length, snp gives an index of a snp, r gives corresponding LD r correlation
    def get_ld_r_tag(self, tag_index):
        num_ld_r = self._check_error(self.cdll.bgmg_num_ld_r_tag(self._context_id, tag_index))
        snp_array = np.zeros(shape=(num_ld_r,), dtype=np.int32)
        r_array = np.zeros(shape=(num_ld_r,), dtype=np.float32)
        self._check_error(self.cdll.bgmg_retrieve_ld_r_tag(self._context_id, tag_index, num_ld_r, snp_array, r_array))
        return (snp_array, r_array)

    # return (snp, tag, r) tuple, representing LD structure of a given chromosome
    # snp, tag and r are vectors of the same length
    # snp gives an index of a reference snp
    # tag gives an index of a tag snp
    # r gives corresponding LD r correlation
    def get_ld_r_chr(self, chr_label):
        num_ld_r = self._check_error(self.cdll.bgmg_num_ld_r_chr(self._context_id, chr_label))
        tag_array = np.zeros(shape=(num_ld_r,), dtype=np.int32)
        snp_array = np.zeros(shape=(num_ld_r,), dtype=np.int32)
        r_array = np.zeros(shape=(num_ld_r,), dtype=np.float32)
        self._check_error(self.cdll.bgmg_retrieve_ld_r_chr(self._context_id, chr_label, num_ld_r, tag_array, snp_array, r_array))
        return (tag_array, snp_array, r_array)

    # return (tag, snp, r) tuple, representing LD structure of a given range of snps
    # [from_tag, to_tag) - left snp inclusive, right snp exclusive, values between 0 and num_snp.
    # tag, snp and r are vectors of the same length
    # tag gives an index of a tag snp
    # snp gives an index of a reference snp
    # r gives corresponding LD r correlation
    def get_ld_r_tag_range(self, from_tag, to_tag):
        num_ld_r = self._check_error(self.cdll.bgmg_num_ld_r_snp_range(self._context_id, from_tag, to_tag))
        snp_array = np.zeros(shape=(num_ld_r,), dtype=np.int32)
        tag_array = np.zeros(shape=(num_ld_r,), dtype=np.int32)
        r_array = np.zeros(shape=(num_ld_r,), dtype=np.float32)
        self._check_error(self.cdll.bgmg_retrieve_ld_r_snp_range(self._context_id, from_tag, to_tag, num_ld_r, tag_array, snp_array, r_array))
        return (tag_array, snp_array, r_array)

    def calc_unified_univariate_aux(self, trait, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL):
        num_component = pi_vec.shape[1]
        num_snp = pi_vec.shape[0]
        aux = np.zeros(shape=(self.num_tag,), dtype=np.float32)
        self.cdll.bgmg_calc_unified_univariate_cost(self._context_id, trait, num_component, num_snp, pi_vec.flatten(), sig2_vec.flatten(), sig2_zeroA, sig2_zeroC, sig2_zeroL, aux)
        self._check_error()
        return aux

    def calc_unified_univariate_cost_with_gradients(self, trait, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL):
        num_component = pi_vec.shape[1]
        num_snp = pi_vec.shape[0]
        aux = np.zeros(shape=(num_snp*num_component+2,), dtype=np.float32)
        cost = self.cdll.bgmg_calc_unified_univariate_cost(self._context_id, trait, num_component, num_snp, pi_vec.flatten(), sig2_vec.flatten(), sig2_zeroA, sig2_zeroC, sig2_zeroL, aux)
        self._check_error()
        return cost, aux

    def calc_unified_univariate_cost(self, trait, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL):
        num_component = pi_vec.shape[1]
        num_snp = pi_vec.shape[0]
        cost = self.cdll.bgmg_calc_unified_univariate_cost(self._context_id, trait, num_component, num_snp, pi_vec.flatten(), sig2_vec.flatten(), sig2_zeroA, sig2_zeroC, sig2_zeroL, None)
        self._check_error()
        return cost

    def calc_unified_univariate_pdf(self, trait, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, zgrid):
        num_component = pi_vec.shape[1]
        num_snp = pi_vec.shape[0]
        zgrid_data = (zgrid if isinstance(zgrid, np.ndarray) else np.array(zgrid)).astype(np.float32)
        pdf = np.zeros(shape=(np.size(zgrid),), dtype=np.float32)
        self._check_error(self.cdll.bgmg_calc_unified_univariate_pdf(self._context_id, trait, num_component, num_snp, pi_vec.flatten(), sig2_vec.flatten(), sig2_zeroA, sig2_zeroC, sig2_zeroL, np.size(zgrid), zgrid_data, pdf))
        return pdf

    def calc_unified_univariate_power(self, trait, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, zthresh, ngrid):
        num_component = pi_vec.shape[1]
        num_snp = pi_vec.shape[0]
        ngrid_data = (ngrid if isinstance(ngrid, np.ndarray) else np.array(ngrid)).astype(np.float32)
        svec_num = np.zeros(shape=(np.size(ngrid),), dtype=np.float32)
        svec_denom = np.zeros(shape=(np.size(ngrid),), dtype=np.float32)
        self._check_error(self.cdll.bgmg_calc_unified_univariate_power(self._context_id, trait, num_component, num_snp, pi_vec.flatten(), sig2_vec.flatten(), sig2_zeroA, sig2_zeroC, sig2_zeroL, zthresh, np.size(ngrid), ngrid_data, svec_num, svec_denom))
        return svec_num, svec_denom

    def calc_unified_univariate_delta_posterior(self, trait, pi_vec, sig2_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL):
        # beta_posterior = delta_posterior / sqrt(H_j N_j)
        c0 = np.zeros(shape=(self.num_tag,), dtype=np.float32)
        c1 = np.zeros(shape=(self.num_tag,), dtype=np.float32)
        c2 = np.zeros(shape=(self.num_tag,), dtype=np.float32)
        num_component = pi_vec.shape[1]
        num_snp = pi_vec.shape[0]
        self._check_error(self.cdll.bgmg_calc_unified_univariate_delta_posterior(self._context_id, trait, num_component, num_snp, pi_vec.flatten(), sig2_vec.flatten(), sig2_zeroA, sig2_zeroC, sig2_zeroL, self.num_tag, c0, c1, c2))
        return (c0, c1, c2)

    def calc_unified_bivariate_cost(self, pi_vec, sig2_beta, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL):
        pi_vec = (pi_vec if isinstance(pi_vec, np.ndarray) else np.array(pi_vec)).astype(np.float32)
        sig2_beta = (sig2_beta if isinstance(sig2_beta, np.ndarray) else np.array(sig2_beta)).astype(np.float32)
        rho_vec = (rho_vec if isinstance(rho_vec, np.ndarray) else np.array(rho_vec)).astype(np.float32)
        sig2_zeroA = (sig2_zeroA if isinstance(sig2_zeroA, np.ndarray) else np.array(sig2_zeroA)).astype(np.float32)
        sig2_zeroC = (sig2_zeroC if isinstance(sig2_zeroC, np.ndarray) else np.array(sig2_zeroC)).astype(np.float32)
        sig2_zeroL = (sig2_zeroL if isinstance(sig2_zeroL, np.ndarray) else np.array(sig2_zeroL)).astype(np.float32)
        aux = np.zeros(shape=(self.num_tag*3,), dtype=np.float32)
        cost = self.cdll.bgmg_calc_unified_bivariate_cost(self._context_id, self.num_snp, pi_vec.flatten(), sig2_beta.flatten(), rho_vec.flatten(), sig2_zeroA.flatten(), sig2_zeroC.flatten(), sig2_zeroL.flatten(), rho_zeroA, rho_zeroL, aux)
        self._check_error()
        return cost

    def calc_unified_bivariate_aux(self, pi_vec, sig2_beta, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL):
        pi_vec = (pi_vec if isinstance(pi_vec, np.ndarray) else np.array(pi_vec)).astype(np.float32)
        sig2_beta = (sig2_beta if isinstance(sig2_beta, np.ndarray) else np.array(sig2_beta)).astype(np.float32)
        rho_vec = (rho_vec if isinstance(rho_vec, np.ndarray) else np.array(rho_vec)).astype(np.float32)
        sig2_zeroA = (sig2_zeroA if isinstance(sig2_zeroA, np.ndarray) else np.array(sig2_zeroA)).astype(np.float32)
        sig2_zeroC = (sig2_zeroC if isinstance(sig2_zeroC, np.ndarray) else np.array(sig2_zeroC)).astype(np.float32)
        sig2_zeroL = (sig2_zeroL if isinstance(sig2_zeroL, np.ndarray) else np.array(sig2_zeroL)).astype(np.float32)
        aux = np.zeros(shape=(self.num_tag*3,), dtype=np.float32)
        cost = self.cdll.bgmg_calc_unified_bivariate_cost(self._context_id, self.num_snp, pi_vec.flatten(), sig2_beta.flatten(), rho_vec.flatten(), sig2_zeroA.flatten(), sig2_zeroC.flatten(), sig2_zeroL.flatten(), rho_zeroA, rho_zeroL, aux)
        self._check_error()
        return aux.reshape([self.num_tag, -1])

    def calc_unified_bivariate_pdf(self, pi_vec, sig2_beta, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL, zvec1, zvec2):
        pi_vec = (pi_vec if isinstance(pi_vec, np.ndarray) else np.array(pi_vec)).astype(np.float32)
        sig2_beta = (sig2_beta if isinstance(sig2_beta, np.ndarray) else np.array(sig2_beta)).astype(np.float32)
        rho_vec = (rho_vec if isinstance(rho_vec, np.ndarray) else np.array(rho_vec)).astype(np.float32)
        sig2_zeroA = (sig2_zeroA if isinstance(sig2_zeroA, np.ndarray) else np.array(sig2_zeroA)).astype(np.float32)
        sig2_zeroC = (sig2_zeroC if isinstance(sig2_zeroC, np.ndarray) else np.array(sig2_zeroC)).astype(np.float32)
        sig2_zeroL = (sig2_zeroL if isinstance(sig2_zeroL, np.ndarray) else np.array(sig2_zeroL)).astype(np.float32)
        zvec1 = (zvec1 if isinstance(zvec1, np.ndarray) else np.array(zvec1)).astype(np.float32)
        zvec2 = (zvec2 if isinstance(zvec2, np.ndarray) else np.array(zvec2)).astype(np.float32)
        if np.size(zvec1) != np.size(zvec2): raise(RuntimeError("len(zvec1) != len(zvec2)"))
        pdf = np.zeros(shape=(np.size(zvec1),), dtype=np.float32)
        self._check_error(self.cdll.bgmg_calc_unified_bivariate_pdf(self._context_id, self.num_snp, pi_vec.flatten(), sig2_beta.flatten(), rho_vec.flatten(), sig2_zeroA.flatten(), sig2_zeroC.flatten(), sig2_zeroL.flatten(), rho_zeroA, rho_zeroL, np.size(zvec1), zvec1, zvec2, pdf))
        return pdf

    def calc_unified_bivariate_delta_posterior(self, pi_vec, sig2_beta, rho_vec, sig2_zeroA, sig2_zeroC, sig2_zeroL, rho_zeroA, rho_zeroL):
        pi_vec = (pi_vec if isinstance(pi_vec, np.ndarray) else np.array(pi_vec)).astype(np.float32)
        sig2_beta = (sig2_beta if isinstance(sig2_beta, np.ndarray) else np.array(sig2_beta)).astype(np.float32)
        rho_vec = (rho_vec if isinstance(rho_vec, np.ndarray) else np.array(rho_vec)).astype(np.float32)
        sig2_zeroA = (sig2_zeroA if isinstance(sig2_zeroA, np.ndarray) else np.array(sig2_zeroA)).astype(np.float32)
        sig2_zeroC = (sig2_zeroC if isinstance(sig2_zeroC, np.ndarray) else np.array(sig2_zeroC)).astype(np.float32)
        sig2_zeroL = (sig2_zeroL if isinstance(sig2_zeroL, np.ndarray) else np.array(sig2_zeroL)).astype(np.float32)
        c00 = np.zeros(shape=(self.num_tag,), dtype=np.float32)
        c10 = np.zeros(shape=(self.num_tag,), dtype=np.float32)
        c01 = np.zeros(shape=(self.num_tag,), dtype=np.float32)
        c20 = np.zeros(shape=(self.num_tag,), dtype=np.float32)
        c11 = np.zeros(shape=(self.num_tag,), dtype=np.float32)
        c02 = np.zeros(shape=(self.num_tag,), dtype=np.float32)
        self._check_error(self.cdll.bgmg_calc_unified_bivariate_delta_posterior(self._context_id, self.num_snp, pi_vec.flatten(), sig2_beta.flatten(), rho_vec.flatten(), sig2_zeroA.flatten(), sig2_zeroC.flatten(), sig2_zeroL.flatten(), rho_zeroA, rho_zeroL, self.num_tag, c00, c10, c01, c20, c11, c02))
        return (c00, c10, c01, c20, c11, c02)

    def annotate_intervals(self, posvec, interval_start, interval_end, interval_label):
        # input:
        #   posvec - coordinates (e.g. base-pair positions) to annotate to a set of intervals
        #   interval_start, interval_end - coordinates where intervals start (inclusive) and end (exclusive)
        #   interval_labels - labels of each interval
        # result:
        #   np.ndarray with shape (nnz, 2), first row gives a list of indices in posvec, second row gives list of labels
        #   this is a sparse matrix of annotations in coo format.
        posvec = (posvec if isinstance(posvec, np.ndarray) else np.array(posvec)).astype(np.int64)
        interval_start = (interval_start if isinstance(interval_start, np.ndarray) else np.array(interval_start)).astype(np.int64)
        interval_end = (interval_end if isinstance(interval_end, np.ndarray) else np.array(interval_end)).astype(np.int64)
        interval_label = (interval_label if isinstance(interval_label, np.ndarray) else np.array(interval_label)).astype(np.int32)

        bytesize = self._check_error(self.cdll.bgmg_request_annotate_intervals(np.size(posvec), posvec.flatten(), np.size(interval_label), interval_start.flatten(), interval_end.flatten(), interval_label.flatten()))
        coo = np.zeros(shape=(bytesize//4,), dtype=np.int32)
        self._check_error(self.cdll.bgmg_copy_requested_object(bytesize, ctypes.c_char_p(coo.ctypes.data)))
        return coo.reshape([2, -1])

    def __str__(self):
        description = []
        for attr_name in '_lib_name', '_context_id', 'num_snp', 'num_tag':
            try:
                attr_value = getattr(self, attr_name)
                description.append('{}: {}'.format(attr_name, attr_value))
            except RuntimeError:
                pass
        return 'LibBgmg({})'.format(', '.join(description))
    __repr__ = __str__

    def _set_vec_impl(self, func, arg_type, arg_value, trait=None):
        data = (arg_value if isinstance(arg_value, np.ndarray) else np.array(arg_value)).astype(arg_type)
        args = [self._context_id] + ([trait] if trait != None else []) + [np.size(data), data]
        self._check_error(func(*args))

    def _get_vec_impl(self, func, arg_type, arg_size, trait=None):
        numpy_ndarray = np.zeros(shape=(arg_size,), dtype=arg_type)
        args = [self._context_id] + ([trait] if trait != None else []) + [np.size(numpy_ndarray), numpy_ndarray]
        self._check_error(func(*args))
        return numpy_ndarray

    def _check_error(self, error_code=None):
        logging.debug('status:{}'.format(self.status))
        
        # if there is an error code, check that it is not negative
        if (error_code is not None) and (error_code < 0):
            raise RuntimeError(self.get_last_error())
        
        # if there is no error code, check that last error message is empty
        if (error_code is None) and (self.get_last_error()):
            raise RuntimeError(self.get_last_error())

        return error_code
  
    def _load_cdll(self, lib_name):
        # choose default library name
        default_lib_name = 'libbgmg.so'
        if sys.platform.startswith('win'):
            default_lib_name = 'bgmg.dll'
        if sys.platform.startswith('darwin'):
            default_lib_name = 'libbgmg.dylib'

        lib_names = []
        
        if lib_name is not None:
            lib_names.append(lib_name)
            
        env_lib_name = os.environ.get('BGMG_SHARED_LIBRARY')
        if env_lib_name is not None:
            lib_names.append(env_lib_name)
        
        lib_names.append(default_lib_name)
        
        # We look into 4 places: lib_name, BGMG_SHARED_LIBRARY, packaged default_lib_name
        # and then default_lib_name
        cdll = None
        exception_message = ""
        for ln in lib_names:
            if not os.path.isfile(ln):
                continue
            try:
                cdll = ctypes.CDLL(ln)
                break
            except OSError as e:
                exception_message += f'{e}'
                continue
        if cdll is None:
            exception_message = (
                '{exception_message}\n'
                'Failed to load BGMG shared library from `{lib_names}`. '
                'Try to add the location of `{default_lib_name}` file into your PATH '
                'system variable, or to set BGMG_SHARED_LIBRARY - the specific system variable '
                'which may point to `{default_lib_name}` file, including the full path.'
            ).format(**locals())
            raise OSError(exception_message)

        return (cdll, ln)

class LibBgmgAnnot(LibBgmg):
    def __init__(self, *args, **kwargs):
        super(LibBgmgAnnot, self).__init__(*args, **kwargs)
        self.annomat = None
        self.annonames = None
        self.genemat = None
        self.genenames = None
        self.bim = None
        self._mafvec_cache = None
        self._tldvec_cache = None

    @property
    def _mafvec(self):
        if self._mafvec_cache is None:
            self._mafvec_cache = self.mafvec
        return self._mafvec_cache

    @property
    def _tldvec(self):
        if self._tldvec_cache is None:
            self._tldvec_cache = self.ld_sum_r2
        return self._tldvec_cache
