# call py.test from <PROJECT_ROOT> folder.
import os, sys, subprocess
sys.path.append(os.getcwd())

import pandas as pd
import numpy as np
import pytest

# features:
# E - extract
# C - complete tag indices
# A - awk a subset of SNPs
tests_overview = '''
    test_ld[21,22]                        - generate --ld-file      (g1000_eur_hm3_chrNN.ld)
    test_snps                             - generate --extract file (g1000_eur_hm3_chr@.snps)
    test_savelib_file                     - generate --loadlib-file (g1000_eur_hm3_chr21to22.bin)
E   test_fit1[1,2]
    test_test1[1,2]
    test_figures_fit1
    test_figures_test1[1,2]
E   test_fit2
    test_test2
    test_figures_fit2

    test_split_sumstats                   - generate --trait1-file  (traitNN.chrNN.sumstats.gz)
    test_plsa_savelib_file                - generate --loadlib-file (g1000_eur_hm3_chrNN.bin)
EC  test_plsa_baseline_trait[1,2]         - generate plsa_traitNN_baseline.[json,weights]
EC  test_plsa_gsa_full_model
EC  test_plsa_model_vs_baseline_enrichment
EC  test_plsa_fit1[1,2]
 CA test_plsa_test1
    test_figures_plsa_fit1
    test_figures_plsa_test1
EC  test_plsa_fit2
 CA test_plsa_test2
    test_figures_plsa_fit2
EC  test_plsa_fit1_infinitesimal[1,2]
EC  test_plsa_fit2_infinitesimal
    test_save_ldsc_reference
'''

data = 'precimed/mixer-test/data'
def subprocess_run(call):
    return subprocess.run(filter(None, ' '.join(call.split('\\')).replace('\n', '').split(' ')))

# py.test precimed/mixer-test/test_cli.py  -k version
def test_version():
    call=f'python precimed/mixer.py --version'
    out = subprocess_run(call)
    assert out.returncode == 0

# py.test precimed/mixer-test/test_cli.py  -k test_ld
@pytest.mark.parametrize("chr_label", [21, 22])
def test_ld(chr_label):
    call=f'python precimed/mixer.py ld --bfile {data}/g1000_eur_hm3_chr{chr_label} --r2min 0.05 --ldscore-r2min 0.01 --lib src/build/lib/libbgmg.so  --out {data}/g1000_eur_hm3_chr{chr_label}.ld --ld-window-kb 10000'
    out = subprocess_run(call)
    assert out.returncode == 0

# py.test precimed/mixer-test/test_cli.py  -k test_snps
def test_snps():
    call=f'python precimed/mixer.py snps --bim-file {data}/g1000_eur_hm3_chr@.bim --ld-file {data}/g1000_eur_hm3_chr@.ld --chr2use 21-22 --subset 20000 --lib src/build/lib/libbgmg.so  --out {data}/g1000_eur_hm3_chr@.snps --seed 123'
    out = subprocess_run(call)
    assert out.returncode == 0

# py.test precimed/mixer-test/test_cli.py  -k test_savelib_file
def test_savelib_file():
    out_file = f'{data}/g1000_eur_hm3_chr21to22.bin'
    if os.path.exists(out_file): os.remove(out_file)
    call=f"""python precimed/mixer.py fit1 \
--lib src/build/lib/libbgmg.so \
--use-complete-tag-indices --exclude-ranges [] \
--savelib-file {data}/g1000_eur_hm3_chr21to22.bin \
--bim-file {data}/g1000_eur_hm3_chr@.bim \
--ld-file {data}/g1000_eur_hm3_chr@.ld \
--chr2use 21-22 \
--out {data}/g1000_eur_hm3_chr21to22
"""
    out = subprocess_run(call)
    assert out.returncode == 0
    assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_fit1
@pytest.mark.parametrize("trait_index", [1, 2])
def test_fit1(trait_index):
    out_file = f'{data}/trait{trait_index}.fit1.json'
    if os.path.exists(out_file): os.remove(out_file)
    call=f"""python precimed/mixer.py fit1 \
--lib src/build/lib/libbgmg.so \
--trait1-file {data}/trait{trait_index}.sumstats.gz \
--bim-file {data}/g1000_eur_hm3_chr@.bim \
--ld-file {data}/g1000_eur_hm3_chr@.ld \
--extract {data}/g1000_eur_hm3_chr@.snps \
--fit-sequence diffevo-fast neldermead-fast \
--chr2use 21-22 --seed 123 --diffevo-fast-repeats 2 \
--cubature-max-evals 5 \
--ci-alpha 0.95 --ci-samples 10 --ci-power-samples 5 \
--out {data}/trait{trait_index}.fit1
"""
    out = subprocess_run(call)
    assert out.returncode == 0
    assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_test1
@pytest.mark.parametrize("trait_index", [1, 2])
def test_test1(trait_index):
    out_file = f'{data}/trait{trait_index}.test1.json'
    if os.path.exists(out_file): os.remove(out_file)
    call=f"""python precimed/mixer.py test1 \
--lib src/build/lib/libbgmg.so \
--trait1-file {data}/trait{trait_index}.sumstats.gz \
--load-params-file {data}/trait{trait_index}.fit1.json \
--bim-file {data}/g1000_eur_hm3_chr@.bim \
--ld-file {data}/g1000_eur_hm3_chr@.ld \
--chr2use 21-22 --seed 123 \
--make-snps-file \
--out {data}/trait{trait_index}.test1
"""
    out = subprocess_run(call)
    assert out.returncode == 0
    assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_figures_fit1
def test_figures_fit1():
    out_file = f'{data}/fit1.csv'
    if os.path.exists(out_file): os.remove(out_file)
    call=f"""python precimed/mixer_figures.py one \
--json {data}/trait1.fit1.json {data}/trait2.fit1.json \
--out {data}/fit1 --trait1 trait_one trait_two --ext svg
"""
    out = subprocess_run(call)
    assert out.returncode == 0    
    assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_figures_test1
@pytest.mark.parametrize("trait_index", [1, 2])
def test_figures_test1(trait_index):
    out_files = [f'{data}/test{trait_index}.csv', f'{data}/test{trait_index}.power.svg', f'{data}/test{trait_index}.power.csv', f'{data}/test{trait_index}.qq.svg', f'{data}/test{trait_index}.qqbin.svg']
    for out_file in out_files:
        if os.path.exists(out_file): os.remove(out_file)
    call=f"""python precimed/mixer_figures.py one \
--json {data}/trait{trait_index}.test1.json \
--out {data}/test{trait_index} --trait1 trait --ext svg
"""
    out = subprocess_run(call)
    assert out.returncode == 0    
    for out_file in out_files: assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_fit2
def test_fit2():
    out_file = f'{data}/fit2.json'
    if os.path.exists(out_file): os.remove(out_file)
    call=f"""python precimed/mixer.py fit2 \
--lib src/build/lib/libbgmg.so \
--trait1-file {data}/trait1.sumstats.gz \
--trait2-file {data}/trait2.sumstats.gz \
--trait1-params-file {data}/trait1.fit1.json \
--trait2-params-file {data}/trait2.fit1.json \
--bim-file {data}/g1000_eur_hm3_chr@.bim \
--ld-file {data}/g1000_eur_hm3_chr@.ld \
--extract {data}/g1000_eur_hm3_chr@.snps \
--chr2use 21-22 --seed 123 --diffevo-fast-repeats 2 \
--fit-sequence diffevo-fast neldermead-fast brute1-fast brent1-fast \
--exclude-ranges chr21:20-21MB chr22:19100-19900KB \
--save-weights \
--out {data}/fit2
"""
    out = subprocess_run(call)
    assert out.returncode == 0
    assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_fit2_region
def test_fit2_region():
    out_file = f'{data}/fit2_region.json'
    if os.path.exists(out_file): os.remove(out_file)
    call=f"""python precimed/mixer.py fit2 \
--lib src/build/lib/libbgmg.so \
--trait1-file {data}/trait1.sumstats.gz \
--trait2-file {data}/trait2.sumstats.gz \
--load-params-file {data}/fit2.json \
--bim-file {data}/g1000_eur_hm3_chr@.bim \
--ld-file {data}/g1000_eur_hm3_chr@.ld \
--extract {data}/g1000_eur_hm3_chr@.snps \
--exclude-ranges chr21:0-18MB chr21:19-300MB chr22:0-25MB chr22:26-300MB \
--chr2use 21-22 --seed 123 \
--fit-sequence neldermead-pi-fast \
--out {data}/fit2_region
"""
    out = subprocess_run(call)
    assert out.returncode == 0
    assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_test2
def test_test2():
    out_file = f'{data}/test2.json'
    if os.path.exists(out_file): os.remove(out_file)
    call=f"""python precimed/mixer.py test2 \
--lib src/build/lib/libbgmg.so \
--trait1-file {data}/trait1.sumstats.gz \
--trait2-file {data}/trait2.sumstats.gz \
--load-params-file {data}/fit2.json \
--bim-file {data}/g1000_eur_hm3_chr@.bim \
--ld-file {data}/g1000_eur_hm3_chr@.ld \
--chr2use 21-22 --seed 123 \
--hardprune-r2 0.1 \
--downsample-factor 1000 \
--kmax-pdf 5 --kmax 50 \
--make-snps-file \
--out {data}/test2
"""
    out = subprocess_run(call)
    assert out.returncode == 0
    assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_figures_fit2
def test_figures_fit2():
    out_files = [f'{data}/fit2.csv', f'{data}/fit2.svg']
    for out_file in out_files:
        if os.path.exists(out_file): os.remove(out_file)
    
    call=f"""python precimed/mixer_figures.py two \
--json-fit {data}/fit2.json --json-test {data}/test2.json \
--out {data}/fit2 --trait1 trait_one --trait2 trait_two --ext svg
"""
    out = subprocess_run(call)
    assert out.returncode == 0    
    for out_file in out_files: assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_split_sumstats
def test_split_sumstats():
    for trait in ['trait1', 'trait2']:
        out_files = [f'{data}/{trait}.chr21.sumstats.gz', f'{data}/{trait}.chr22.sumstats.gz']
        for out_file in out_files:
            if os.path.exists(out_file): os.remove(out_file)
        call = f"python precimed/mixer.py split_sumstats --trait1-file precimed/mixer-test/data/{trait}.sumstats.gz --out precimed/mixer-test/data/{trait}.chr@.sumstats.gz --chr2use 21-22"
        out = subprocess_run(call)
        assert out.returncode == 0
        for out_file in out_files: 
            assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_plsa_savelib_file
def test_plsa_savelib_file():
    out_files = [f'{data}/g1000_eur_hm3_chr{chri}.bin' for chri in [21,22]]
    for out_file in out_files:
        if os.path.exists(out_file): os.remove(out_file)
    call=f"""python precimed/mixer.py plsa \
--lib src/build/lib/libbgmg.so \
--use-complete-tag-indices --exclude-ranges [] \
--savelib-file {data}/g1000_eur_hm3_chr@.bin \
--bim-file {data}/g1000_eur_hm3_chr@.bim \
--ld-file {data}/g1000_eur_hm3_chr@.ld \
--chr2use 21-22 \
--out {data}/g1000_eur_hm3_chr@
"""
    out = subprocess_run(call)
    assert out.returncode == 0
    for out_file in out_files:
        assert os.path.exists(out_file)

def run_plsa(baseline, trait, se_samples=None, go_all_genes_label=None):
    out_files = ['.json', '.done']
    if not baseline: out_files += ['.go_test_enrich.csv']
    out_files = [f"{data}/{f'plsa_{trait}_baseline' if baseline else f'plsa_{trait}_model'}{file}" for file in out_files]
    for out_file in out_files:
        if os.path.exists(out_file): os.remove(out_file)

    call=f"""python precimed/mixer.py plsa \
--lib src/build/lib/libbgmg.so \
--trait1-file {data}/{trait}.chr@.sumstats.gz \
--use-complete-tag-indices \
--bim-file {data}/g1000_eur_hm3_chr@.bim \
--loadlib-file {data}/g1000_eur_hm3_chr@.bin \
--annot-file {data}/g1000_eur_hm3_chr@.annot.gz \
--go-file {data}/{"go-file-baseline" if baseline else "go-file-gene"}.csv \
--extract {data}/g1000_eur_hm3_chr@.snps \
--exclude-ranges chr21:20-21MB chr22:19100-19900KB \
--chr2use 21-22 --seed 123 \
--adam-epoch 3 3 --adam-step 0.064 0.032 \
--out {data}/{f"plsa_{trait}_baseline" if baseline else f"plsa_{trait}_model"}"""
    call += ' --gsa-base' if baseline else ' --gsa-full'
    call += '' if se_samples is None else f' --se-samples {se_samples} '
    call += '' if go_all_genes_label is None else f' --go-all-genes-label {go_all_genes_label} '
    if not baseline: call += f" --load-params-file {data}/plsa_trait1_baseline.json"
    if not baseline: call += f" --go-file-test {data}/go-file-geneset.csv "
    out = subprocess_run(call)
    assert out.returncode == 0
    for out_file in out_files: 
        assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_plsa_baseline_trait1
def test_plsa_baseline_trait1():
    run_plsa(baseline=True, trait='trait1')

# py.test precimed/mixer-test/test_cli.py  -k test_plsa_baseline_trait2
def test_plsa_baseline_trait2():
    run_plsa(baseline=True, trait='trait2')

# py.test precimed/mixer-test/test_cli.py  -k test_plsa_gsa_full_model_withSE
def test_plsa_gsa_full_model_withSE():
    run_plsa(baseline=False, trait='trait1')

# py.test precimed/mixer-test/test_cli.py  -k test_plsa_gsa_full_model_woSE
def test_plsa_gsa_full_model_woSE():
    run_plsa(baseline=False, trait='trait1', se_samples=0)

# py.test precimed/mixer-test/test_cli.py  -k test_plsa_gsa_full_model_justbase
def test_plsa_gsa_full_model_justbase():
    run_plsa(baseline=False, trait='trait1', se_samples=None, go_all_genes_label='no_such_geneset')

# py.test precimed/mixer-test/test_cli.py  -k test_plsa_model_vs_baseline_enrichment
def run_plsa_model_vs_baseline_enrichment(calc_loglike_diff_go_test):
    trait = 'trait1'
    out_files = ['.json', '.done', '.annot_test_enrich.csv', '.go_test_enrich.csv']
    out_files = [f"{data}/plsa_{trait}_model-vs-baseline{file}" for file in out_files]
    for out_file in out_files:
        if os.path.exists(out_file): os.remove(out_file)

    call=f"""python precimed/mixer.py plsa \
--lib src/build/lib/libbgmg.so \
--trait1-file {data}/{trait}.chr@.sumstats.gz \
--use-complete-tag-indices \
--bim-file {data}/g1000_eur_hm3_chr@.bim \
--loadlib-file {data}/g1000_eur_hm3_chr@.bin \
--annot-file {data}/g1000_eur_hm3_chr@.annot.gz \
--annot-file-test {data}/g1000_eur_hm3_chr@.annot.gz \
--go-file {data}/go-file-gene.csv \
--go-file-test {data}/go-file-geneset.csv \
--adam-disable \
--load-params-file {data}/plsa_trait1_model.json \
--load-baseline-params-file {data}/plsa_trait1_baseline.json \
--calc-loglike-diff-go-test {calc_loglike_diff_go_test} \
--extract {data}/g1000_eur_hm3_chr@.snps \
--exclude-ranges chr21:20-21MB chr22:19100-19900KB \
--chr2use 21-22 --seed 123 \
--out {data}/plsa_{trait}_model-vs-baseline"""
    out = subprocess_run(call)
    assert out.returncode == 0
    for out_file in out_files:
        assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_plsa_model_vs_baseline_enrichment_fast_calc
def test_plsa_model_vs_baseline_enrichment_fast_calc():
    run_plsa_model_vs_baseline_enrichment(calc_loglike_diff_go_test='fast')

# py.test precimed/mixer-test/test_cli.py  -k test_plsa_model_vs_baseline_enrichment_full_calc
def test_plsa_model_vs_baseline_enrichment_full_calc():
    run_plsa_model_vs_baseline_enrichment(calc_loglike_diff_go_test='full')

# py.test precimed/mixer-test/test_cli.py  -k test_plsa_fit1
@pytest.mark.parametrize("trait_index", [1, 2])
def test_plsa_fit1(trait_index):
    out_file = f'{data}/plsa_trait{trait_index}_baseline.fit1.json'
    if os.path.exists(out_file): os.remove(out_file)
    call=f"""python precimed/mixer.py fit1 \
--lib src/build/lib/libbgmg.so \
--trait1-file {data}/trait{trait_index}.sumstats.gz \
--load-params-file {data}/plsa_trait{trait_index}_baseline.json \
--use-complete-tag-indices \
--bim-file {data}/g1000_eur_hm3_chr@.bim \
--loadlib-file {data}/g1000_eur_hm3_chr21to22.bin \
--annot-file {data}/g1000_eur_hm3_chr@.annot.gz \
--go-file {data}/go-file-baseline.csv \
--extract {data}/g1000_eur_hm3_chr@.snps \
--exclude-ranges chr21:20-21MB chr22:19100-19900KB \
--chr2use 21-22 --seed 123 --diffevo-fast-repeats 2 \
--cubature-max-evals 5 \
--fit-sequence diffevo-fast neldermead-fast \
--ci-alpha 0.95 --ci-samples 10 \
--out {data}/plsa_trait{trait_index}_baseline.fit1
"""
    out = subprocess_run(call)
    assert out.returncode == 0
    assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_plsa_test1
def test_plsa_test1():
    out_file = f'{data}/plsa_trait1_baseline.test1.json'
    if os.path.exists(out_file): os.remove(out_file)
    os.system(f"zcat {data}/trait1.sumstats.gz | awk 'NR % 2 == 1' | gzip > {data}/trait1.filtered.sumstats.gz")
    call=f"""python precimed/mixer.py test1 \
--lib src/build/lib/libbgmg.so \
--trait1-file {data}/trait1.filtered.sumstats.gz \
--load-params-file {data}/plsa_trait1_baseline.fit1.json \
--use-complete-tag-indices \
--bim-file {data}/g1000_eur_hm3_chr@.bim \
--ld-file {data}/g1000_eur_hm3_chr@.ld \
--annot-file {data}/g1000_eur_hm3_chr@.annot.gz \
--go-file {data}/go-file-baseline.csv \
--chr2use 21-22 --seed 123 \
--make-snps-file --kmax 1000 \
--out {data}/plsa_trait1_baseline.test1
"""
    out = subprocess_run(call)
    assert out.returncode == 0
    assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_figures_plsa_fit1
def test_figures_plsa_fit1():
    out_file = f'{data}/plsa_trait1_baseline.fit1.csv'
    if os.path.exists(out_file): os.remove(out_file)
    call=f"""python precimed/mixer_figures.py one \
--json {data}/plsa_trait1_baseline.fit1.json  \
--out {data}/plsa_trait1_baseline.fit1 --trait1 trait_one --ext svg
"""
    out = subprocess_run(call)
    assert out.returncode == 0    
    assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_figures_plsa_test1
def test_figures_plsa_test1():
    out_files = [f'{data}/plsa_trait1_baseline.test1.csv',
                 f'{data}/plsa_trait1_baseline.test1.power.svg',
                 f'{data}/plsa_trait1_baseline.test1.power.csv',
                 f'{data}/plsa_trait1_baseline.test1.qq.svg',
                 f'{data}/plsa_trait1_baseline.test1.qqbin.svg']
    for out_file in out_files:
        if os.path.exists(out_file): os.remove(out_file)
    call=f"""python precimed/mixer_figures.py one \
--json {data}/trait1.test1.json \
--out {data}/plsa_trait1_baseline.test1 --trait1 trait_one --ext svg
"""
    out = subprocess_run(call)
    assert out.returncode == 0    
    for out_file in out_files: assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_plsa_fit2
def test_plsa_fit2():
    out_file = f'{data}/plsa_baseline.fit2.json'
    if os.path.exists(out_file): os.remove(out_file)
    call=f"""python precimed/mixer.py fit2 \
--lib src/build/lib/libbgmg.so \
--trait1-file {data}/trait1.sumstats.gz \
--trait2-file {data}/trait2.sumstats.gz \
--trait1-params-file {data}/plsa_trait1_baseline.fit1.json \
--trait2-params-file {data}/plsa_trait2_baseline.fit1.json \
--use-complete-tag-indices \
--bim-file {data}/g1000_eur_hm3_chr@.bim \
--ld-file {data}/g1000_eur_hm3_chr@.ld \
--annot-file {data}/g1000_eur_hm3_chr@.annot.gz \
--go-file {data}/go-file-baseline.csv \
--extract {data}/g1000_eur_hm3_chr@.snps \
--chr2use 21-22 --seed 123 --diffevo-fast-repeats 2 \
--fit-sequence diffevo-fast neldermead-fast brute1-fast brent1-fast \
--exclude-ranges chr21:20-21MB chr22:19100-19900KB \
--out {data}/plsa_baseline.fit2
"""
    out = subprocess_run(call)
    assert out.returncode == 0
    assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_plsa_test2
def test_plsa_test2():
    out_file = f'{data}/plsa_baseline.test2.json'
    if os.path.exists(out_file): os.remove(out_file)
    os.system(f"zcat {data}/trait1.sumstats.gz | awk 'NR % 2 == 1' | gzip > {data}/trait1.filtered.sumstats.gz")
    os.system(f"zcat {data}/trait2.sumstats.gz | awk 'NR % 2 == 1' | gzip > {data}/trait2.filtered.sumstats.gz")
    call=f"""python precimed/mixer.py test2 \
--lib src/build/lib/libbgmg.so \
--trait1-file {data}/trait1.filtered.sumstats.gz \
--trait2-file {data}/trait2.filtered.sumstats.gz \
--load-params-file {data}/plsa_baseline.fit2.json \
--use-complete-tag-indices \
--bim-file {data}/g1000_eur_hm3_chr@.bim \
--ld-file {data}/g1000_eur_hm3_chr@.ld \
--annot-file {data}/g1000_eur_hm3_chr@.annot.gz \
--go-file {data}/go-file-baseline.csv \
--chr2use 21-22 --seed 123 \
--downsample-factor 100 --kmax-pdf 10 \
--out {data}/plsa_baseline.test2
"""
    out = subprocess_run(call)
    assert out.returncode == 0
    assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_figures_plsa_fit2
def test_figures_plsa_fit2():
    out_files = [f'{data}/plsa_baseline_fit2.csv', f'{data}/plsa_baseline_fit2.svg']
    for out_file in out_files:
        if os.path.exists(out_file): os.remove(out_file)
    
    call=f"""python precimed/mixer_figures.py two \
--json-fit {data}/plsa_baseline.fit2.json --json-test {data}/plsa_baseline.test2.json \
--out {data}/plsa_baseline_fit2 --trait1 trait_one --trait2 trait_two --ext svg
"""
    out = subprocess_run(call)
    assert out.returncode == 0    
    for out_file in out_files: assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_plsa_fit1_infinitesimal
@pytest.mark.parametrize("trait_index", [1,2])
def test_plsa_fit1_infinitesimal(trait_index):
    out_file = f'{data}/plsa_trait{trait_index}_infinitesimal.fit1.json'
    if os.path.exists(out_file): os.remove(out_file)
    call=f"""python precimed/mixer.py fit1 \
--lib src/build/lib/libbgmg.so \
--trait1-file {data}/trait{trait_index}.sumstats.gz \
--load-params-file {data}/plsa_trait{trait_index}_baseline.json \
--use-complete-tag-indices \
--bim-file {data}/g1000_eur_hm3_chr@.bim \
--ld-file {data}/g1000_eur_hm3_chr@.ld \
--annot-file {data}/g1000_eur_hm3_chr@.annot.gz \
--go-file {data}/go-file-baseline.csv \
--extract {data}/g1000_eur_hm3_chr@.snps \
--chr2use 21-22 --seed 123 --diffevo-fast-repeats 2 \
--cubature-max-evals 5 \
--fit-sequence infinitesimal \
--out {data}/plsa_trait{trait_index}_infinitesimal.fit1
"""
    out = subprocess_run(call)
    assert out.returncode == 0
    assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_plsa_fit2_infinitesimal
def test_plsa_fit2_infinitesimal():
    out_file = f'{data}/plsa_infinitesimal.fit2.json'
    if os.path.exists(out_file): os.remove(out_file)
    call=f"""python precimed/mixer.py fit2 \
--lib src/build/lib/libbgmg.so \
--trait1-file {data}/trait1.sumstats.gz \
--trait2-file {data}/trait2.sumstats.gz \
--trait1-params-file {data}/plsa_trait1_baseline.json \
--trait2-params-file {data}/plsa_trait2_baseline.json \
--use-complete-tag-indices \
--bim-file {data}/g1000_eur_hm3_chr@.bim \
--ld-file {data}/g1000_eur_hm3_chr@.ld \
--annot-file {data}/g1000_eur_hm3_chr@.annot.gz \
--go-file {data}/go-file-baseline.csv \
--extract {data}/g1000_eur_hm3_chr@.snps \
--chr2use 21-22 --seed 123 --diffevo-fast-repeats 2 \
--fit-sequence infinitesimal \
--out {data}/plsa_infinitesimal.fit2
"""
    out = subprocess_run(call)
    assert out.returncode == 0
    assert os.path.exists(out_file)

# py.test precimed/mixer-test/test_cli.py  -k test_save_ldsc_reference
def test_save_ldsc_reference():
    call=f"""python precimed/mixer.py save_ldsc_reference \
--lib src/build/lib/libbgmg.so \
--use-complete-tag-indices \
--bim-file {data}/g1000_eur_hm3_chr@.bim \
--annot-file {data}/g1000_eur_hm3_chr@.annot.gz \
--loadlib-file {data}/g1000_eur_hm3_chr@.bin \
--annot-file {data}/g1000_eur_hm3_chr@.annot.gz \
--chr2use 21-22 --seed 123 \
--save-ldsc-reference \
--out {data}/ldsc_reference_chr@
"""
    out = subprocess_run(call)
    assert out.returncode == 0
    assert os.path.exists(f'{data}/ldsc_reference_chr21.mat')
    assert os.path.exists(f'{data}/ldsc_reference_chr22.mat')
