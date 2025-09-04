import numpy as np
import numdifftools as nd
import scipy
from scipy.sparse import coo_matrix
from scipy.sparse import csr_matrix

from common.libbgmg import _auxoption_ezvec2, _auxoption_gradients, _auxoption_none, _auxoption_tagpdf, _auxoption_tagpdferr

from common.utils import _log_exp_converter
from common.utils import _logit_logistic_converter
from common.utils import make_one_zeros_vec

# params for MAF-, LD-, and annotation-dependent architectures
# note that AnnotUnivariateParams.sig2_beta parameter is different from UnivariateParams.sig2_beta, more specifically:
# UnivariateParams.sig2_beta is equivalent to AnnotUnivariateParams.sig2_beta / AnnotUnivariateParams.pi
# This means that UnivariateParametrization_natural_axis works more or less like (now obsolete) UnivariateParametrization
class AnnotUnivariateParams(object):
    def __init__(self, pi=None, sig2_beta=None, sig2_annot=None, sig2_gene=None, s=None, l=None, sig2_zeroA=None, sig2_zeroL=None, libbgmg=None, snp_defvec=None, p_annot=1, p_gene=1):
        self._pi = pi
        self._sig2_beta = sig2_beta
        self._sig2_annot = np.array(sig2_annot).flatten() if (sig2_annot is not None) else make_one_zeros_vec(len(libbgmg.annonames))
        self._sig2_gene = np.array(sig2_gene).flatten() if (sig2_gene is not None) else make_one_zeros_vec(len(libbgmg.genenames))
        self._s = s
        self._l = l
        self._sig2_zeroA = sig2_zeroA
        self._sig2_zeroL = sig2_zeroL
        self._libbgmg = libbgmg
        self._snp_defvec = snp_defvec
        self._p_annot = p_annot  # define p in (x_1^p + ... + x_t^p)^(1/p) for aggregating sigma^2_t annotations functional annotation and gene sets
        self._p_gene = p_gene    # p=1 means sum; p->inf is similar to max; in practice values between p=2 and p=5 are good choices.

    def copy(self, pi=None, sig2_beta=None, sig2_zeroA=None, libbgmg=None):
        return AnnotUnivariateParams(
            pi=self._pi if (pi is None) else pi,
            sig2_beta=self._sig2_beta if (sig2_beta is None) else sig2_beta,
            sig2_annot=self._sig2_annot,
            sig2_gene=self._sig2_gene,
            s=self._s,
            l=self._l,
            sig2_zeroA=self._sig2_zeroA if (sig2_zeroA is None) else sig2_zeroA,
            sig2_zeroL=self._sig2_zeroL,
            libbgmg=self._libbgmg if (libbgmg is None) else libbgmg,
            snp_defvec=self._snp_defvec,
            p_annot=self._p_annot,
            p_gene=self._p_gene)

    @property
    def pi(self):
        return self._pi

    @property
    def sig2_beta(self):
        return self._sig2_beta

    @property
    def sig2_zeroA(self):
        return self._sig2_zeroA

    @property
    def _mafvec(self):
        return self._libbgmg._mafvec

    @property
    def _tldvec(self):
        return self._libbgmg._tldvec

    @property
    def _annomat(self):
        return self._libbgmg.annomat

    @property
    def _annonames(self):
        return self._libbgmg.annonames

    @property
    def _genemat(self):
        return self._libbgmg.genemat

    @property
    def _genenames(self):
        return self._libbgmg.genenames

    @property
    def _n_annot(self):
        return len(self._sig2_annot)

    @property
    def _n_genes(self):
        return len(self._sig2_gene)

    @property
    def _genenames(self):
        return self._libbgmg.genenames

    def as_string(self, attrs_list=['_pi', '_sig2_beta', '_s', '_l', '_sig2_zeroA', '_sig2_zeroL', '_p_annot', '_p_gene', "_n_annot", "_n_genes"]):
        description = []
        for attr_name in attrs_list:
            try:
                attr_value = getattr(self, attr_name)
                description.append('{}: {}'.format(attr_name, attr_value))
            except RuntimeError:
                pass
        return 'AnnotUnivariateParams({})'.format(', '.join(description))

    def __str__(self):
        return self.as_string()

    __repr__ = __str__

    def find_pi_mat(self, num_snp):
        return np.matmul(np.ones(shape=(num_snp, 1), dtype=np.float32), np.array(self._pi, dtype=np.float32).reshape(1, -1))

    def find_sig2_mat_without_gene(self, sig2_annot=None):
        hetvec = np.float32(2.0) * self._mafvec * (1-self._mafvec)
        het_pow_s = np.zeros(hetvec.shape, dtype=np.float32)
        het_pow_s[hetvec > 0] = np.power(hetvec[hetvec > 0], np.float32(self._s))
        het_pow_s = het_pow_s / np.mean(het_pow_s)

        tld_pow_l = np.zeros(self._tldvec.shape, dtype=np.float32)
        tld_pow_l[self._tldvec > 0] = np.power(self._tldvec[self._tldvec > 0], np.float32(self._l))
        tld_pow_l = tld_pow_l / np.mean(tld_pow_l)

        sig2_beta_effective = self.sig2_beta / self.pi  # divide by pi, to ensure that pi has no effect on heritability (as this helps to improve convergence)
        #print('sig2_beta_effective: ', sig2_beta_effective)

        sig2_annot_array = np.maximum(0, np.array(self._sig2_annot if (sig2_annot is None) else sig2_annot)).flatten()
        sig2_annot_scale = np.mean(sig2_annot_array)
        sig2_annot_p = np.power(sig2_annot_array / sig2_annot_scale, self._p_annot).astype(np.float32)

        sig2_mat = np.maximum(0, np.matmul(
                np.multiply(
                    sig2_annot_scale * np.power(self._annomat.dot(sig2_annot_p), 1/self._p_annot),
                    np.multiply(het_pow_s, tld_pow_l)).reshape(-1, 1),
                np.array(sig2_beta_effective, dtype=np.float32).reshape(1, -1)))
        if self._snp_defvec is not None: sig2_mat = np.multiply(sig2_mat, self._snp_defvec)

        if not np.all(np.isfinite(sig2_mat)): raise(ValueError('Not finite sig2_mat'))
        #print('sig2mat[:10] = ', sig2_mat[:10])
        return sig2_mat

    def find_sig2_mat(self, sig2_annot=None, sig2_gene=None):
        # shape:
        # sig2_beta (ncat, )
        # _annomat (nsnp, nannot)
        # _genemat (nsnp, ngene)
        # _sig2_annot (nannot, )
        # _sig2_gene (ngene, )
        # _mafvec (nsnp, )
        # _tldvec (nsnp, )
        # return: (nsnp, ncat)

        hetvec = np.float32(2.0) * self._mafvec * (1-self._mafvec)
        het_pow_s = np.zeros(hetvec.shape, dtype=np.float32)
        het_pow_s[hetvec > 0] = np.power(hetvec[hetvec > 0], np.float32(self._s))
        het_pow_s = het_pow_s / np.mean(het_pow_s)

        tld_pow_l = np.zeros(self._tldvec.shape, dtype=np.float32)
        tld_pow_l[self._tldvec > 0] = np.power(self._tldvec[self._tldvec > 0], np.float32(self._l))
        tld_pow_l = tld_pow_l / np.mean(tld_pow_l)

        sig2_beta_effective = self.sig2_beta / self.pi  # divide by pi, to ensure that pi has no effect on heritability (as this helps to improve convergence)
        #print('sig2_beta_effective: ', sig2_beta_effective)

        sig2_annot_array = np.maximum(0, np.array(self._sig2_annot if (sig2_annot is None) else sig2_annot, dtype=np.float32)).flatten()
        sig2_annot_scale = np.mean(sig2_annot_array)
        sig2_annot_p = np.power(sig2_annot_array / sig2_annot_scale, self._p_annot).astype(np.float32)
        
        sig2_gene_array = np.maximum(0, np.array(self._sig2_gene if (sig2_gene is None) else sig2_gene, dtype=np.float32)).flatten()
        sig2_gene_scale = np.mean(sig2_gene_array)
        sig2_gene_p = np.power(sig2_gene_array / sig2_gene_scale, self._p_gene).astype(np.float32)

        sig2_mat = np.maximum(0, np.matmul(
                np.multiply(
                    np.multiply(sig2_annot_scale * np.power(self._annomat.dot(sig2_annot_p), 1/self._p_annot),
                                sig2_gene_scale * np.power(self._genemat.dot(sig2_gene_p), 1/self._p_gene)),
                            np.multiply(het_pow_s, tld_pow_l)).reshape(-1, 1),
                np.array(sig2_beta_effective, dtype=np.float32).reshape(1, -1)))
        if self._snp_defvec is not None: sig2_mat = np.multiply(sig2_mat, self._snp_defvec)

        if not np.all(np.isfinite(sig2_mat)): raise(ValueError('Not finite sig2_mat'))
        #print('sig2mat[:10] = ', sig2_mat[:10])
        return sig2_mat

    def find_h2(self):
        h2 = self.find_annot_h2(annomat=np.ones(shape=(self._libbgmg.num_snp, 1), dtype=np.float32))
        return h2.flat[0]

    def find_h2_vec(self):
        sig2_mat = self.find_sig2_mat()
        hetvec = np.float32(2.0) * self._mafvec * (1-self._mafvec)
        h2_vec = np.multiply(hetvec, np.matmul(sig2_mat, self.find_pi_mat(num_snp=1).reshape([-1, 1])).flatten())
        return h2_vec

    def find_annot_h2(self, annomat=None):
        h2_vec = self.find_h2_vec()
        return (annomat if (annomat is not None) else self._annomat).transpose().dot(h2_vec.reshape((len(h2_vec), 1))).reshape([1, -1])

    # find the proportion of causal SNPs expected to explain a given percentage of SNP-based heritability
    def find_nckoef(self, quantile=0.9, num_iter=10000):
        if self._pi == 1.0: return np.nan
        sig2_mat = self.find_sig2_mat().flatten()
        hetvec = np.float32(2.0) * self._mafvec * (1-self._mafvec)

        rng = np.random.default_rng()
        result = np.empty((num_iter, ))
        for i in range(num_iter):
            num_samples = np.random.binomial(n=len(hetvec), p=self._pi)
            if num_samples==0:
                result[i] = 0
                continue
            idx = rng.choice(len(sig2_mat), num_samples, replace=False, shuffle=False)
            beta = np.multiply(np.random.standard_normal(num_samples), np.sqrt(sig2_mat[idx]) )
            h2_vec = np.multiply(np.power(beta, 2), hetvec[idx])
            h2_vec = np.sort(h2_vec)[::-1]
            h2_vec = np.cumsum(h2_vec)
            result[i] = np.argmax(h2_vec >= (quantile * h2_vec[-1]))
        NCKoef = np.mean(result) / (len(hetvec) * self._pi)
        NCKoef_se = (np.std(result) / (len(hetvec) * self._pi)) / np.sqrt(len(result))
        self._libbgmg.log_message(f"Extimated NCKoef={NCKoef} (SE={NCKoef_se})")
        return NCKoef

    def find_mean_sig2_beta(self):
        sig2_mat = self.find_sig2_mat().flatten()
        return np.mean(sig2_mat)

    def find_annot_snps(self, annomat=None):
        return np.sum(annomat if (annomat is not None) else self._annomat, 0)

    def find_annot_enrich(self, annomat=None):
        h2_annot = self.find_annot_h2(annomat)
        snps_annot = self.find_annot_snps(annomat)
        h2_total = h2_annot[0][0]
        snps_total = snps_annot[0]
        return np.divide(np.divide(h2_annot, h2_total), np.divide(snps_annot, snps_total))

    def cost(self, lib, trait):
        pi_mat = self.find_pi_mat(lib.num_snp)
        sig2_mat = self.find_sig2_mat()
        lib.set_option('aux_option', _auxoption_none)
        value = lib.calc_unified_univariate_cost(trait, pi_mat, sig2_mat, self._sig2_zeroA, sig2_zeroC=1, sig2_zeroL=self._sig2_zeroL)
        return value if np.isfinite(value) else 1e100

    def cost_with_gradients(self, lib, trait):
        pi_mat = self.find_pi_mat(lib.num_snp)
        sig2_mat = self.find_sig2_mat()
        lib.set_option('aux_option', _auxoption_gradients)
        value, aux = lib.calc_unified_univariate_cost_with_gradients(trait, pi_mat, sig2_mat, self._sig2_zeroA, sig2_zeroC=1, sig2_zeroL=self._sig2_zeroL)
        if self._snp_defvec is not None: aux[:-2] = np.multiply(aux[:-2], np.squeeze(self._snp_defvec))
        return value if np.isfinite(value) else 1e100, aux

    def pdf(self, lib, trait, zgrid):
        pi_mat = self.find_pi_mat(lib.num_snp)
        sig2_mat = self.find_sig2_mat()
        return lib.calc_unified_univariate_pdf(trait, pi_mat, sig2_mat, self._sig2_zeroA, sig2_zeroC=1, sig2_zeroL=self._sig2_zeroL, zgrid=zgrid)

    def tag_pdf(self, lib, trait):
        pi_mat = self.find_pi_mat(lib.num_snp)
        sig2_mat = self.find_sig2_mat()
        lib.set_option('aux_option', _auxoption_tagpdf)
        retval = lib.calc_unified_univariate_aux(trait, pi_mat, sig2_mat, self._sig2_zeroA, sig2_zeroC=1, sig2_zeroL=self._sig2_zeroL)
        lib.set_option('aux_option', _auxoption_none)
        return retval

    def tag_pdf_err(self, lib, trait):
        pi_mat = self.find_pi_mat(lib.num_snp)
        sig2_mat = self.find_sig2_mat()
        lib.set_option('aux_option', _auxoption_tagpdferr)
        retval = lib.calc_unified_univariate_aux(trait, pi_mat, sig2_mat, self._sig2_zeroA, sig2_zeroC=1, sig2_zeroL=self._sig2_zeroL)
        lib.set_option('aux_option', _auxoption_none)
        return retval

    def tag_ez2(self, lib, trait):
        pi_mat = self.find_pi_mat(lib.num_snp)
        sig2_mat = self.find_sig2_mat()
        lib.set_option('aux_option', _auxoption_ezvec2)
        retval = lib.calc_unified_univariate_aux(trait, pi_mat, sig2_mat, self._sig2_zeroA, sig2_zeroC=1, sig2_zeroL=self._sig2_zeroL)
        lib.set_option('aux_option', _auxoption_none)
        return retval

    def power(self, lib, trait, ngrid, zthresh=5.45):
        pi_mat = self.find_pi_mat(lib.num_snp)
        sig2_mat = self.find_sig2_mat()
        svec_num, svec_denom = lib.calc_unified_univariate_power(trait, pi_mat, sig2_mat, self._sig2_zeroA, sig2_zeroC=1, sig2_zeroL=self._sig2_zeroL, zthresh=zthresh, ngrid=ngrid)
        return np.divide(svec_num, svec_denom)

    def delta_posterior(self, lib, trait):
        pi_mat = self.find_pi_mat(lib.num_snp)
        sig2_mat = self.find_sig2_mat()
        return lib.calc_unified_univariate_delta_posterior(trait, pi_mat, sig2_mat, self._sig2_zeroA, sig2_zeroC=1, sig2_zeroL=self._sig2_zeroL)

class AnnotUnivariateParametrization(object):
    def __init__(self, trait, constraint, sig2_gene_base=None):
        self._trait = trait
        self._constraint = constraint # of type AnnotUnivariateParams, None indicate files that must be searched
        for constraint_s2b in self._constraint._sig2_annot:
            if (constraint_s2b is not None) and (not np.isfinite(constraint_s2b)):
                raise(ValueError('sig2_annot must be None or finite (nan or inf is not accepted)'))
        for constraint_s2b in self._constraint._sig2_gene:
            if (constraint_s2b is not None) and (not np.isfinite(constraint_s2b)):
                raise(ValueError('sig2_gene must be None or finite (nan or inf is not accepted)'))

        # hack-hack, separate trick to handle --constrain-base-gene-category
        self._sig2_gene_base = sig2_gene_base

        # hack-hack, separate trick to handle --nullify-go-all-genes-sig2-gene
        self._nullify_indices = []

        # hack-hack, trick to speed up computations w.r.t. sig2_gene vector
        self._B = None
        self._gene_indices_on_chr = None
        self._enable_faster_gradient_computation_chr = None

    @property
    def _lib(self):
        return self._constraint._libbgmg

    def params_to_vec(self, params):
        vec = []
        for constraint_s2b, params_s2b in zip(self._constraint._sig2_annot, params._sig2_annot):
            if constraint_s2b is None: vec.append(_log_exp_converter(params_s2b, invflag=False))
        for constraint_s2b, params_s2b in zip(self._constraint._sig2_gene, params._sig2_gene):
            if constraint_s2b is None: vec.append(_log_exp_converter(params_s2b, invflag=False))
        for constraint_pi, params_pi in zip([self._constraint._pi], [params._pi]):
            if constraint_pi is None: vec.append(_logit_logistic_converter(params_pi, invflag=False))
        for constraint_s2b, params_s2b in zip([self._constraint._sig2_beta], [params._sig2_beta]):
            if constraint_s2b is None: vec.append(_log_exp_converter(params_s2b, invflag=False))
        if self._constraint._s is None: vec.append(params._s)
        if self._constraint._l is None: vec.append(params._l)
        if self._constraint._sig2_zeroA is None: vec.append(_log_exp_converter(params._sig2_zeroA, invflag=False))
        if self._constraint._sig2_zeroL is None: vec.append(_log_exp_converter(params._sig2_zeroL, invflag=False))
        return [float(x) for x in vec]

    def vec_to_params(self, vec):
        vec = list(vec)
        sig2_annot_value = [(constraint_s2b if (constraint_s2b is not None) else _log_exp_converter(vec.pop(0), invflag=True)) for constraint_s2b in self._constraint._sig2_annot]
        sig2_gene_value = [(constraint_s2b if (constraint_s2b is not None) else _log_exp_converter(vec.pop(0), invflag=True)) for constraint_s2b in self._constraint._sig2_gene]
        pi_vec = [(constraint_pi if (constraint_pi is not None) else _logit_logistic_converter(vec.pop(0), invflag=True)) for constraint_pi in [self._constraint._pi]]
        sig2_beta_vec = [(constraint_s2b if (constraint_s2b is not None) else _log_exp_converter(vec.pop(0), invflag=True)) for constraint_s2b in [self._constraint._sig2_beta]]

        if self._sig2_gene_base is not None:
            sig2_gene_value[0] = self._sig2_gene_base
        for idx in self._nullify_indices:
            sig2_gene_value[idx] = 0

        return AnnotUnivariateParams(
            pi=pi_vec[0],
            sig2_beta=sig2_beta_vec[0],
            sig2_annot=sig2_annot_value,
            sig2_gene=sig2_gene_value,
            s=self._constraint._s if (self._constraint._s is not None) else vec.pop(0),
            l=self._constraint._l if (self._constraint._l is not None) else vec.pop(0),
            sig2_zeroA=self._constraint._sig2_zeroA if (self._constraint._sig2_zeroA is not None) else _log_exp_converter(vec.pop(0), invflag=True),
            sig2_zeroL=self._constraint._sig2_zeroL if (self._constraint._sig2_zeroL is not None) else _log_exp_converter(vec.pop(0), invflag=True),
            libbgmg=self._constraint._libbgmg,
            p_annot=self._constraint._p_annot,
            p_gene=self._constraint._p_gene)

    def find_sig2_jacobian(self, sig2_mat, sig2_annot, p_annot, annomat):
        # sig2_annot_p                            - col-vector, Tx1, T=number of annotations
        # sig2_mat and sig2_annot_snps_factor     - col-vector, Mx1, M=number of SNPs
        
        # for numerical accuracy we re-scale sig2_annot (so it's mean is 1.0) - this way we can still use np.power(sig2_annot, p_annot), even with np.float32
        sig2_annot_array = np.maximum(0, np.array(sig2_annot)).flatten()
        sig2_annot_scale = np.mean(sig2_annot_array)
        sig2_annot_p = np.power(sig2_annot_array / sig2_annot_scale, p_annot).astype(np.float32).reshape(-1, 1)

        with np.errstate(invalid='ignore', divide='ignore'):
            sig2_annot_snps_factor = np.multiply(sig2_mat, np.power(annomat.dot(sig2_annot_p), -1).reshape(-1, 1))
            sig2_annot_snps_factor[~np.isfinite(sig2_annot_snps_factor)] = 0

        M = len(sig2_annot_snps_factor)
        T = len(sig2_annot_p)

        sig2_annot_p_matrix = coo_matrix((sig2_annot_p.flat, (np.arange(T), np.arange(T))), shape=(T,T), dtype=np.float32).tocsr()
        sig2_annot_snps_matrix = coo_matrix((sig2_annot_snps_factor.flat, (np.arange(M), np.arange(M))), shape=(M,M), dtype=np.float32).tocsr()

        sig2_annot_jacob = (sig2_annot_snps_matrix.dot(annomat)).dot(sig2_annot_p_matrix)
        return sig2_annot_jacob

    def find_sig2_annot_jacobian(self, params):
        sig2_annot_jacob = None; sig2_gene_jacob = None

        if (None not in self._constraint._sig2_annot) and (None not in self._constraint._sig2_gene):
            return sig2_annot_jacob, sig2_gene_jacob

        sig2_mat = params.find_sig2_mat()

        if None in self._constraint._sig2_annot:
            sig2_annot_jacob = self.find_sig2_jacobian(sig2_mat, params._sig2_annot, params._p_annot, params._annomat)

        if None in self._constraint._sig2_gene:
            sig2_gene_jacob = self.find_sig2_jacobian(sig2_mat, params._sig2_gene, params._p_gene, params._genemat)

        return sig2_annot_jacob, sig2_gene_jacob

    def find_jacobian(self, params):
        # find jacobian matrix (jacob[i, j]) of partial derivatives where
        #   - index 'i' runs across elements in vec (those that we optimize)
        #   - index 'j' runs across elements in aux vector, calculated from 'aux_option=4', with length (num_components*num_snps + 2)
        # This require constrained pi. Also, the code below also relies on assumption that link function between optimization space and
        # sig2_beta is an exponent, f(x)=e^x, which somewhat simplifies expressions for partial derivatives 

        sig2_mat = params.find_sig2_mat()
        num_components = sig2_mat.shape[1]
        aux_len = sig2_mat.size + 2

        jacob = []

        for constraint_pi, _ in zip([self._constraint._pi], [params._pi]):
            if constraint_pi is None:
                grad=np.zeros((aux_len, ))
                grad[:] = np.nan
                jacob.append(grad)

        for constraint_s2b, _ in zip([self._constraint._sig2_beta], [params._sig2_beta]):
            if constraint_s2b is None:
                grad=np.zeros((aux_len, ))
                grad[0:-2] = sig2_mat.flatten()
                jacob.append(grad)

        if self._constraint._s is None:
            hetvec = np.float32(2.0) * params._mafvec * (1-params._mafvec)
            log_hetvec = np.zeros(hetvec.shape, dtype=np.float32)
            log_hetvec[hetvec > 0] = np.log(hetvec[hetvec > 0])

            het_pow_s = np.zeros(hetvec.shape, dtype=np.float32)
            het_pow_s[hetvec > 0] = np.power(hetvec[hetvec > 0], np.float32(params._s))
            het_pow_s = het_pow_s / np.mean(het_pow_s)

            total_het_pow_s = np.sum(het_pow_s)
            total_het_pow_s_deriv = np.sum(np.multiply(log_hetvec, het_pow_s))

            log_hetvec_with_offset = log_hetvec - total_het_pow_s_deriv / total_het_pow_s

            grad=np.zeros((aux_len, ))
            grad[0:-2] = np.multiply(sig2_mat, log_hetvec_with_offset.reshape([-1, 1]) * np.ones((1, num_components))).flatten()
            jacob.append(grad)

        if self._constraint._l is None:
            log_tldvec = np.zeros(params._tldvec.shape, dtype=np.float32)
            log_tldvec[params._tldvec > 0] = np.log(params._tldvec[params._tldvec > 0])

            tld_pow_l = np.zeros(params._tldvec.shape, dtype=np.float32)
            tld_pow_l[params._tldvec > 0] = np.power(params._tldvec[params._tldvec > 0], np.float32(params._l))
            tld_pow_l = tld_pow_l / np.mean(tld_pow_l)

            total_tld_pow_l = np.sum(tld_pow_l)
            total_tld_pow_l_deriv = np.sum(np.multiply(log_tldvec, tld_pow_l))

            log_tldvec_with_offset = log_tldvec - total_tld_pow_l_deriv / total_tld_pow_l

            grad=np.zeros((aux_len, ))
            grad[0:-2] = np.multiply(sig2_mat, log_tldvec_with_offset.reshape([-1, 1]) * np.ones((1, num_components))).flatten()
            jacob.append(grad)

        if self._constraint._sig2_zeroA is None:
            grad=np.zeros((aux_len, ))
            grad[-2] = params._sig2_zeroA
            jacob.append(grad)

        if self._constraint._sig2_zeroL is None:
            grad=np.zeros((aux_len, ))
            grad[-1] = params._sig2_zeroL
            jacob.append(grad)

        return np.array(jacob).astype(np.float32)

    def calc_cost(self, vec, verbose=True):
        params = self.vec_to_params(vec)
        if verbose: self._lib.log_message(params.as_string())
        return params.cost(self._lib, self._trait)

    def calc_faster_gradient_computation(libbgmg, params, chr_label, go_all_genes_label):
        # enable faster gradients computation w.r.g. sig2gene parameters
        # under simplifying conditions: pi==1, z1max=Inf, p_a = g_a = 1 (additive overlap model)
        tag_array, snp_array, r_array = libbgmg.get_ld_r_chr(chr_label)
        r2=csr_matrix((np.power(r_array, 2), (tag_array, snp_array)), dtype=np.float32)
        hetvec = scipy.sparse.diags(2 * np.multiply(libbgmg._mafvec, 1-libbgmg._mafvec), dtype=np.float32)
        nvec1 = libbgmg.nvec1
        nvec1[~np.isfinite(nvec1)] = 0
        nvec1 = scipy.sparse.diags(nvec1, dtype=np.float32)

        A2 = (nvec1.dot(r2)).dot(hetvec)

        gene_indices_on_chr = np.where(libbgmg.genemat.sum(0) > 0)[1]
        skip_first_rows = 0
        if libbgmg.genenames[gene_indices_on_chr[0]] in ['base', go_all_genes_label]:
            skip_first_rows = 1
            if libbgmg.genenames[gene_indices_on_chr[1]] in ['base', go_all_genes_label]:
                assert libbgmg.genenames[gene_indices_on_chr[0]] == 'base'
                assert libbgmg.genenames[gene_indices_on_chr[1]] == go_all_genes_label
                skip_first_rows = 2
        gene_indices_on_chr = gene_indices_on_chr[skip_first_rows:]

        S = scipy.sparse.diags(params.find_sig2_mat_without_gene().flatten(), dtype=np.float32)
        G = S.dot(libbgmg.genemat[:, gene_indices_on_chr])
        B = A2.dot(G)

        return B, gene_indices_on_chr

    def enable_faster_gradient_computation(self, params, libbgmg, chr_label, go_all_genes_label):
        # this option is intended for scenario where we don't need gradient for "base" and "coding_genes" categories
        assert self._sig2_gene_base is not None
        assert len(self._nullify_indices) == 1
        params._libbgmg = libbgmg
        self._constraint._libbgmg = libbgmg
        self._B, self._gene_indices_on_chr = AnnotUnivariateParametrization.calc_faster_gradient_computation(
            libbgmg, params, chr_label, go_all_genes_label)
        self._enable_faster_gradient_computation_chr = chr_label

    def disable_faster_gradient_computation(self):
        self._B = None
        self._gene_indices_on_chr = None
        self._enable_faster_gradient_computation_chr = None

    def calc_cost_with_gradients(self, vec, vec_step, verbose=True, force_numeric=False):
        params = self.vec_to_params(vec)
        if verbose: self._lib.log_message(params.as_string())

        gradients_fast = None
        if (not force_numeric) and (self._gene_indices_on_chr is not None):
            sig2_gene_chr = params._sig2_gene[self._gene_indices_on_chr]
            sig2_vec = self._B.dot(sig2_gene_chr) + params.sig2_zeroA

            #tag_pdf = 0.5 * np.log(2*np.pi * sig2_vec) + 0.5 * np.divide(np.power(self._lib.zvec1, 2), sig2_vec)
            #tag_pdf = tag_pdf[self._lib.weights>0]
            #print('tag_pdf_fast:', tag_pdf, tag_pdf.shape)

            cost_fast = np.nansum(
                np.multiply(self._lib.weights,
                        0.5 * np.log(2*np.pi * sig2_vec) + 0.5 * np.divide(np.power(self._lib.zvec1, 2), sig2_vec)))

            weights_gradient = np.multiply(self._lib.weights, np.divide(np.power(self._lib.zvec1, 2) - sig2_vec, 2 * np.power(sig2_vec, 2)))
            gradient_sig2_vec = self._B.transpose().dot(weights_gradient)  # loglike gradient w.r.t. sig2_gene;
            gradient_xvec = np.multiply(gradient_sig2_vec, sig2_gene_chr)  # loglike gradient w.r.t. xvec, (where sig2_gene = e^xvec).
            gradients_fast = np.zeros((len(vec), 1))
            gradients_fast[self._gene_indices_on_chr] = gradient_xvec.reshape([-1, 1])

            gradients = gradients_fast
            cost = cost_fast

        elif not force_numeric:
            #tag_pdf = -np.log(params.tag_pdf(self._lib, self._trait))
            #tag_pdf = tag_pdf[self._lib.weights>0]
            #print('tag_pdf     :', tag_pdf, tag_pdf.shape)

            cost, aux = params.cost_with_gradients(self._lib, self._trait)
            #print(cost, ', '.join([str(x) for x in aux]))
            if not np.all(np.isfinite(aux)): raise(ValueError('nan value present in aux, can not compute gradients'))
            sig2_annot_jacobian, sig2_gene_jacobian = self.find_sig2_annot_jacobian(params)
            jacob = self.find_jacobian(params)
            sig2_annot_gradients = None if (sig2_annot_jacobian is None) else sig2_annot_jacobian.transpose().dot(aux[:-2].reshape([-1, 1]))
            sig2_gene_gradients = None if (sig2_gene_jacobian is None) else sig2_gene_jacobian.transpose().dot(aux[:-2].reshape([-1, 1]))
            gradients = np.matmul(jacob, aux.reshape([-1, 1])) if jacob.size else None
            gradients = np.concatenate([x for x in [sig2_annot_gradients, sig2_gene_gradients, gradients] if (x is not None)])
            #if gradients_fast is not None:
            #    print('cost_fast: ', cost_fast, 'gradients_fast:', ', '.join([str(x) for x in gradients_fast[2:]]))
            #    print('cost     : ', cost,      'gradients     :', ', '.join([str(x) for x in gradients[2:]]))
            #    print('\n')

        else:
            cost = params.cost(self._lib, self._trait)
            gradients = np.zeros((len(vec), 1)); gradients[:] = np.nan

        for index, gradient in enumerate(gradients):
            if np.isfinite(gradient): continue
            if vec_step is None: raise ValueError('Unable to compute numeric gradients, vec_step is None')
            vec_plus = np.array(vec); vec_plus[index] = vec_plus[index] + vec_step[index]
            vec_minus = np.array(vec); vec_minus[index] = vec_minus[index] - vec_step[index]
            cost_plus = self.calc_cost(vec_plus)
            cost_minus = self.calc_cost(vec_minus)
            gradients[index] = -1 * (cost_plus - cost_minus) / (2 * vec_step[index])  # flipped sign (going towards minimum, not maximum)
        return cost, gradients

def _params_to_dict(params):
    if isinstance(params, AnnotUnivariateParams):
        return {'pi': params.pi, 'sig2_beta': params.sig2_beta,
                'sig2_zeroA': params._sig2_zeroA, 'sig2_zeroL': params._sig2_zeroL,
                's': params._s, 'l': params._l,
                'sig2_annot': params._sig2_annot, 'annonames': params._annonames, 'p_annot': params._p_annot,
                'sig2_gene': params._sig2_gene, 'genenames': params._genenames, 'p_gene': params._p_gene,
                'type':'AnnotUnivariateParams' }
    else:
        raise ValueError('Unrecognized type: {}'.format(type(params)))

def find_lib_sig2_beta(lib_catnames, json_catnames, json_sig2_cats, default=0):
    # check duplicates, filter out zeros
    assert(len(lib_catnames) == len(set(lib_catnames)))
    assert(len(json_catnames) == len(set(json_catnames)))
    json_catnames, json_sig2_cats = zip(*[(n, b) for (n, b) in zip(json_catnames, json_sig2_cats) if (b != 0)])

    # check if some of json cats (with non-zero sig2_beta) aren't present in --annot-file / --go-file
    missing = list(set(json_catnames).difference(set(lib_catnames)))
    if len(missing) > 0: raise(ValueError('Categories "' + ", ".join(missing) + '" from --load-params-file are missing in --annot-file / --go-file'))

    mapping = dict(zip(json_catnames, json_sig2_cats))
    return [mapping.get(k, default) for k in lib_catnames]

def _dict_to_params(p, libbgmg, args=None):
    if 'type' not in p:
        raise ValueError('Unrecognized type in dict_to_params()')
    if p['type'] == 'AnnotUnivariateParams':
        sig2_annot = find_lib_sig2_beta(libbgmg.annonames, p['annonames'], p['sig2_annot']) if libbgmg else p['sig2_annot']

        default = 0
        if (args is not None) and ('go_all_genes_label' in args):
            genenames_to_sig2_gene = dict(zip(p['genenames'], p['sig2_gene']))
            if (args.go_all_genes_label in genenames_to_sig2_gene):
                default = genenames_to_sig2_gene[args.go_all_genes_label]
                libbgmg.log_message(f'default sig2_beta={default} ("{args.go_all_genes_label}" category) is used')
            elif 'base' in genenames_to_sig2_gene:
                # this logic makes some sense because 'base' is defined as a completent to all other gene categories
                default = genenames_to_sig2_gene['base']
                libbgmg.log_message(f'default sig2_beta={default} ("base" category) is used')
            else:
                raise ValueError(f'"base" gene-set category is missing in --load-params-file/--load-baseline-params-file')

        sig2_gene = find_lib_sig2_beta(libbgmg.genenames, p['genenames'], p['sig2_gene'], default) if libbgmg else p['sig2_gene']
        return AnnotUnivariateParams(
            pi=p['pi'], s=p['s'], l=p['l'], sig2_beta=p['sig2_beta'],
            sig2_zeroA=p['sig2_zeroA'], sig2_zeroL=p['sig2_zeroL'],
            sig2_annot=sig2_annot, p_annot=p['p_annot'] if ('p_annot' in p) else 1.0,
            sig2_gene=sig2_gene, p_gene=p['p_gene'] if ('p_gene' in p) else 1.0,
            libbgmg=libbgmg)
    else:
        raise ValueError('Unrecognized type: {}'.format(p['type']))
