'''
module load CMake/3.15.3-GCCcore-8.3.0
module load Boost/1.73.0-GCCcore-8.3.0
module load Python/3.7.4-GCCcore-8.3.0
source /cluster/projects/p697/users/ofrei/py3/bin/activate
'''

import pandas as pd
import sys

_base_complement = {"A":"T", "C":"G", "G":"C", "T":"A"}

alleles = ['A', 'T', 'C', 'G']
#matched_pairs = 'ACAC, ACTG, ACCA, ACGT, AGAG, AGTC, AGCT, AGGA, TCAG, TCTC, TCCT, TCGA, TGAC, TGTG, TGCA, TGGT, CAAC, CATG, CACA, CAGT, CTAG, CTTC, CTCT, CTGA, GAAG, GATC, GACT, GAGA, GTAC, GTTG, GTCA, GTGT'.split(', ')
matched_pairs = 'ACAC, ACCA, AGAG, AGGA, TCTC, TCCT, TGTG, TGGT, CAAC, CACA, CTTC, CTCT, GAAG, GAGA, GTTG, GTGT'.split(', ')

if __name__ == "__main__":
    chri = sys.argv[1]  # 21
    prefix = sys.argv[2] # PGC_SCZ_0518_EUR.chr21
    ss = pd.read_csv(prefix, sep='\t')
    ss = ss[ss['A1'].isin(alleles) & ss['A2'].isin(alleles)].copy()
    bim = pd.read_csv(f'/ess/p697/cluster/projects/moba_qc_imputation/resources/HRC/plink/hrc_chr{chri}_EUR_qc_keep1k.bim', delim_whitespace=True, header=None, names='CHR SNP GP BP A1 A2'.split())
    bim = bim[bim['A1'].isin(alleles) & bim['A2'].isin(alleles)].copy()
    df = pd.merge(ss[['SNP', 'A1', 'A2']], bim[['SNP', 'A1', 'A2']], on='SNP', how='inner')
    print(f'{len(df)} SNPs overlap' )
    df['alleles'] = df['A1_x'] + df['A2_x'] + df['A1_y'] + df['A2_y']
    df = df[df['alleles'].isin(matched_pairs)]
    df[['SNP']].to_csv(f'{prefix}.overlapHRC.justrs', sep='\t', index=False, header=None)
    print(f'{len(df)} SNPs match' )
