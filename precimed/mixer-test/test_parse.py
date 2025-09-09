# call py.test from <PROJECT_ROOT> folder.
import os, sys
sys.path.append(os.getcwd())

import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix
from precimed.common import libbgmg
import random

data = 'precimed/mixer-test/data'

_base_complement = {"A":"T", "C":"G", "G":"C", "T":"A"}
def _complement(seq):
    return "".join([_base_complement[b] for b in seq])

def _reverse_complement_variant(variant):
    # variant should be a 2-elemet sequence with upper case string elements
    return ("".join([_base_complement[b] for b in variant[0][::-1]]),
            "".join([_base_complement[b] for b in variant[1][::-1]]))

# py.test precimed/mixer-test/test_parse.py  -k test_parse
def test_parse():
    np.random.seed(123); random.seed(123)
    chr2use = [21, 22]
    bim_file = f"{data}/g1000_eur_hm3_chr@.bim"
    trait1_file = f"{data}/random.sumstats.gz"

    letters = ["A", "T", "C", "G"]
    for chri in chr2use:
        bim = pd.read_csv(f"{data}/g1000_eur_hm3_chr{chri}.bim", sep='\t', header=None, names='CHR SNP GP BP A1 A2'.split())
        if True:  # generate multiallelic (for all SNPs)
            bim['A2'] = ["".join(random.choices(letters, k=random.randint(2, 10))) for _ in range(len(bim)) ]
        if True:  # generate single-allelic code (unless A1==A2, in which case keep original)
            bim['A2_copy'] = ["".join(random.choices(letters, k=1)) for _ in range(len(bim)) ]
            idx = bim['A2_copy'] != bim['A1']
            bim.loc[idx, 'A2'] = bim.loc[idx, 'A2_copy']
            del bim['A2_copy']
        bim.to_csv(f"{data}/g1000_eur_hm3_chr{chri}_multi.bim", sep='\t', index=False, header=None)
        bim_file=f"{data}/g1000_eur_hm3_chr@_multi.bim"

    ref=pd.concat([pd.read_csv(bim_file.replace('@', str(chr_label)), sep='\t', header=None, names='CHR SNP GP BP A1 A2'.split()) for chr_label in chr2use])
    num_snp = len(ref)
    ref['N'] = np.arange(1e4, 1e4+num_snp, 1)  # at this point it's not ref - it's a sumstats
    ref['Z'] = np.random.normal(size=(num_snp,))

    if True:
        # randomly remove half of the SNPs
        idx = np.random.uniform(size=(num_snp,)) < 0.5
        ref.loc[idx, 'Z'] = np.nan

    gtypes = list(zip(ref.A1.apply(str.upper),ref.A2.apply(str.upper)))
    ind_ambiguous = [ gt == _reverse_complement_variant(gt)[::-1] for gt in gtypes]

    true_zscore = ref['Z'].values.copy()
    true_zscore[ind_ambiguous] = np.nan
    #ref.loc[ind_ambiguous, 'Z'] = np.nan  # keep ambiguous Z scores - those should be detected by bgmg_parse (that's what we're testing)

    if True:
        # randomly swap alleles (and flip Z-score)
        idx = np.random.uniform(size=(num_snp,)) < 0.5
        ref.loc[idx, 'A1'], ref.loc[idx, 'A2'], ref.loc[idx, 'Z'] = ref.loc[idx, 'A2'], ref.loc[idx, 'A1'], -ref.loc[idx, 'Z']    

    if True:
        # randomly replace SNPs with complement (and don't flip Z score)
        idx = np.random.uniform(size=(num_snp,)) <= 0.5
        ref.loc[idx, 'A1'] = [_complement(x)[::-1] for x in ref.loc[idx, 'A1'].values]
        ref.loc[idx, 'A2'] = [_complement(x)[::-1] for x in ref.loc[idx, 'A2'].values]

    if True:
        # delete value in SNP column  to force merging on CHR:BP
        idx = np.random.uniform(size=(num_snp,)) <= 0.5
        ref.loc[idx, 'SNP'] = "."

        # introduce non-numeric codes in CHR:BP (which should be ignored if match on SNP is possible)
        ref['BP'] = ref['BP'].astype(str)
        ref.loc[~idx, 'BP'] = 'X'

    if True:
        # reshuffle rows
        ref = ref.sample(frac=1, random_state=42) 

    ref.to_csv(trait1_file, sep='\t', index=False)

    lib=libbgmg.LibBgmg('src/build/lib/libbgmg.so', init_log='precimed/mixer-test/test.log', context_id=0, dispose=True)
    lib.set_option('use_complete_tag_indices', True)
    lib.init(bim_file, "", [21, 22], trait1_file, "", "", "", "")

    zvec1 = lib.zvec1
    #print('\n')
    #print(np.concatenate([zvec1[:10].reshape((1, -1)), true_zscore[:10].reshape((1, -1))]).T)
    #print(ref.head(10))

    idx = np.isfinite(true_zscore)
    assert(np.all(idx == np.isfinite(zvec1)))
    assert(np.all(np.isclose(true_zscore[idx], zvec1[idx])))
