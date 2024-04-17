## Introduction

GSA-MiXeR is a new technique for competitive gene-set analysis, which fits a model for gene-set heritability enrichments for complex human traits, thus allowing the quantification of partitioned heritability and fold enrichment for small gene-sets.

This tutorial has the following parts:

* [Install GSA-MiXeR](#install-gsa-mixer) using Docker or singularity containers;
* [Getting Started Example](#getting-started-example) using tiny dummy data
* [Input Data Formats](#input-data-formats) for summary statistics, functional annotations and gene-set definitions
* [Download LD Reference Files](#download-ld-reference-files) build from for 1kG, UKB or HRC genotypes
* [(optional) Generate LD Reference](#optional-generate-ld-reference) using your own genotype panel
* [Perform GSA-MiXeR and MAGMA analyses](#perform-gsa-mixer-and-magma-analyses) on real-world data
* [Output formats of GSA-MiXeR analysis](#output-formats-of-gsa-mixer-analysis)
* [Command-line reference](#command-line-reference)

## Install GSA-MiXeR

GSA-MiXeR is available as pre-compiled singularity container which can be downloaded from [here](https://github.com/norment/ofrei_repo).

Singularity software (https://sylabs.io/docs/) is most likely available in your HPC environment, however it's also
not too hard to get it up an running on your laptop (especially on Ubuntu, probably also on older MAC with an intel CPU).
To install singularity on Ubuntu follow steps described here: https://sylabs.io/guides/3.7/user-guide/quick_start.html
Note that ``sudo apt-get`` can give only a very old version of singularity, which isn't sufficient.
Therefore it's best to build singularity locally. 
Note that building singularity from source code depends on [GO](https://go.dev/doc/install), 
so it must be installed first. One you have singularity up and running, it might be usefult o have a look at
["singularity shell" options](https://sylabs.io/guides/3.2/user-guide/cli/singularity_shell.html#options) and
[Bind paths and mounts](https://sylabs.io/guides/3.2/user-guide/bind_paths_and_mounts.html) pages from the documentation.

The ``mixer.sif`` container is based on the following [Dockerfile](Dockerfile), which is included in this repository.
In order to do this you will also need to install Docker (e.g. [here](https://docs.docker.com/desktop/install/ubuntu/) are instructions for installing it on Ubuntu). This is not necessary if you are planning to use pre-compiled ``mixer.sif`` container.
The only limitation with pre-compiled ``mixer.sif`` is that it only works for intel-based CPUs, and does not support
ARM architectures such as for example newer Macbook laptops with M1 chip.
But if you would like to build ``mixer.sif`` container from scratch, then the steps are as follows:

```
# check Docker software is installed
>docker --version
Docker version 20.10.7, build 20.10.7-0ubuntu5~21.04.2

# check singularity software is installed
>singularity --version
singularity version 3.7.4

# produce mixer.sif
>docker build -t mixer -f Dockerfile . && scripts/from_docker_image.sh mixer
```

## Getting Started Example

This section depends on example data files located in ``precimed/mixer-test/data`` folder of this repository:
* ```g1000_eur_hm3_chr21to22.[bed,bim,fam]``` - EUR subset of 1kG Phase3 individuals (N=503) for M=34958 SNPs from chr21 and chr22, already constrained to HapMap3 SNPs
*  ```g1000_eur_hm3_chr[21,22].ld``` - LD matrix derived from the above genotypes using ``mixer.py ld`` command
* ```g1000_eur_hm3_chr@.snps``` - SNPs used to subset GWAS z-scores used in fit procedure; the set of SNPs is derived from the above genotypes with ``mixer.py snps`` command
* ```trait1.sumstats.gz``` and ```trait2.sumstats.gz``` - GWAS summary statistics for two traits (only the first trait is used in GSA-MiXeR demo example)
* ```partial.pheno``` - two synthesized phenotypes each with SNP-h2=0.7, generated from the above genotypes using additive genetic model; this file is not used by GSA-MiXeR, but it was used to produce the above GWAS summary statistics via ``plink2 --glm`` call.
* ```g1000_eur_hm3_chr[21,22].annot.gz``` - randomly generated functional annotations in sLDSC format
* ```go-file-baseline.csv``` - baseline model with three gene sets (all_genes, coding_genes, pseudo_genes);
* ```go-file-gene.csv``` - enrichment model with in total 435 real genes from chr21 and chr22
* ```go-file-geneset.csv``` - enrichment model with 562 real gene-sets (constrained to genes on chr21 and chr22)

GSA-MiXeR usage can be illustrated with the following steps.
The first two steps (``$MIXER ld`` and ``$MIXER snps``) are optional, as they relate to producing reference files, which for real-data analysis usually should be downloaded via the links provided below.
Note that ``@`` symbol must remain as it is in all commands, i.e. you don't need to exchange it with a specific chromosome label.
All commands below assume that demo data is locate in your current folder.
Expected execution time of all commands below on a standard laptop is less than 60 seconds.

```
# point MIXER_SIF variable to the location of the mixer.sif file
export MIXER_SIF=<ROOT>/mixer.sif

# define MIXER_PY command which executes mixer.py script from within singularity container
export MIXER_PY="singularity exec --home $PWD:/home ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"

for chri in {21..22}; do ${MIXER_PY} ld --bfile g1000_eur_hm3_chr$chri --r2min 0.05 --ldscore-r2min 0.01 --out g1000_eur_hm3_chr$chri.ld --ld-window-kb 10000; done  

${MIXER_PY} snps --bim-file g1000_eur_hm3_chr@.bim --ld-file g1000_eur_hm3_chr@.ld --chr2use 21-22 --r2 0.6 --maf 0.05 --subset 20000 --out g1000_eur_hm3_chr@.snps --seed 123

# split summary statistics into one file per chromosome
${MIXER_PY} split_sumstats --trait1-file trait1.sumstats.gz --out trait1.chr@.sumstats.gz --chr2use 21-22

# fit baseline model, and use it to calculate heritability attributed to gene-sets in go-file-geneset.csv
${MIXER_PY} plsa --gsa-base \
--trait1-file trait1.chr@.sumstats.gz \
--use-complete-tag-indices \
--bim-file g1000_eur_hm3_chr@.bim \
--ld-file g1000_eur_hm3_chr@.ld \
--annot-file g1000_eur_hm3_chr@.annot.gz \
--go-file go-file-baseline.csv \
--go-file-test go-file-geneset.csv \
--extract g1000_eur_hm3_chr@.snps \
--exclude-ranges chr21:20-21MB chr22:19100-19900KB \
--chr2use 21-22 --seed 123 \
--adam-epoch 3 3 --adam-step 0.064 0.032 \
--out plsa_baseline

# fit enrichment model, and use it to calculate heritability attributed to gene-sets in go-file-geneset.csv
${MIXER_PY} plsa --gsa-full \
--trait1-file trait1.chr@.sumstats.gz \
--use-complete-tag-indices \
--bim-file g1000_eur_hm3_chr@.bim \
--ld-file g1000_eur_hm3_chr@.ld \
--annot-file g1000_eur_hm3_chr@.annot.gz \
--go-file go-file-gene.csv \
--go-file-test go-file-geneset.csv \
--extract g1000_eur_hm3_chr@.snps \
--load-params-file plsa_baseline.json \
--exclude-ranges chr21:20-21MB chr22:19100-19900KB \
--chr2use 21-22 --seed 123 \
--adam-epoch 3 3 --adam-step 0.064 0.032 \
--out plsa_model
```

The commands above are customized to run the analysis faster.
For real-data analysis the commands will need the following changes:
* remove ``--exclude-ranges chr21:20-21MB chr22:19100-19900KB``; by default ``--exclude-ranges`` will exclude MHC region
* remove ``--chr2use 21-22``; by default ``--chr2use`` applies to all chromosomes
* remove ``--adam-epoch 3 3 --adam-step 0.064 0.032``, as this stops Adam fit procedure too early

## Input Data Formats

GSA-MiXeR format for summary statistics (``--trait1-file``) is compatible with LD Score Regression
(i.e. the ``.sumstats.gz`` files), and must include the following columns:
* Eithe one of the following:
  * ``SNP`` or ``RSID`` (marker name), or
  * ``CHR`` (chromosome label) and  ``BP`` or ``POS`` (genomic corrdinates), in a build that is compatible with the reference build (``--bim-file`` argument)
* ``A1`` or ``EffectAllele`` (reference allele)
* ``A2`` or ``OtherAllele`` (alternative allele)
* ``N`` (sample size); for case-control studies this should be the effective sample size computed as ``N=4/(1/ncases+1/ncontrols)``
* ``Z`` (signed test statistic)
Column names must be exactly as defined above, except for upper/lower case which can be arbitrary (all column names from the input file are converted to lower case prior to matching them with expected column names defined above).

It's beneficial to have both ``SNP`` and ``CHR``/``BP`` columns in the data.
In this situation matching SNPs with the reference (``--bim-file``) is performed on marker name (``SNP`` column).
For the remaining SNPs GSA-MiXeR attempts to match using CHR:BP:A1:A2 codes, accounting for situations when (1) ``A1``/``A2`` are swapped between GWAS summary statistics and the reference, (2) alleles are coded on different strands, (3) both of the above. 
The sign of z-score is refersed when needed during this procedure.

We advice against filtering down summary statistics to the set of SNPs included in HapMap3.
Prior to running GSA-MiXeR we advice filtering SNPs with bad imputation quality, if ``INFO`` column is available in summary statistics.
If per-SNP sample size is available, we advice filtering out SNPs with N below half of the median sample size across SNPs.
Other filtering options are built into GSA-MiXeR software, including ``--exclude-ranges`` option to filter out special regions such as MHC, and ``--maf`` and ``--randprune-maf`` to filter out based on minor allele frequency.

## Download LD Reference Files

All reference data described below is based on EUR ancestry, and use ``hg19`` / ``GRCh37`` genomic build.

Reference files derived from 1kG Phase3 EUR population are available for download from [here](https://github.com/comorment/mixer/tree/main/reference/ldsc/1000G_EUR_Phase3_plink).

We also have prepared similar reference files derived from [UKB](https://github.com/precimed/mixer_private_docker/tree/main/reference/ukb_EUR_qc) and [HRC](https://github.com/precimed/mixer_private_docker/tree/main/reference/hrc_EUR_qc), which we plan to release together with GSA-MiXeR tool. Functional annotations are derived from [sLDSC baselineLD_v2.2](https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/baselineLD_v2.2_bedfiles.tgz) using adapted scripts from [here](https://github.com/mkanai/eas_partitioned_ldscore) to annotate our UKB and HRC references. Note that one does not need to compute LD-scores for these annotations, because MiXeR does this internally using sparse LD matrix stored in ``--ld-file`` it receives as an argument.

```
ukb_EUR_qc/about_UKB_qc.txt                                 # overview of QC procedure
ukb_EUR_qc/ukb_imp_chr[1-22]_v3_qc.bim                      # ``--bim-file`` argument
ukb_EUR_qc/baseline_v2.2_ukb_imp_chr[1-22]_v3_qc.annot.gz   # ``--annot-file`` / ``--annot-file-test`` arguments
ukb_EUR_qc/ukb_imp_chr[1-22]_v3_qc.run1.ld                  # ``--ld-file`` argument
ukb_EUR_qc/ukb_imp_chr@_qc.prune_rand2M_rep[1-20].snps      # ``--extract`` argument

hrc_EUR_qc/about_HRC_qc.txt                                 # overview of QC procedure
hrc_EUR_qc/hrc_chr[1-22]_EUR_qc.bim                         # ``--bim-file`` argument
hrc_EUR_qc/hrc_chr[1-22]_EUR_qc.run1.ld                     # ``--annot-file`` / ``--annot-file-test`` arguments
hrc_EUR_qc/baseline_v2.2_hrc_chr[1-22]_EUR_qc.annot.gz      # ``--ld-file`` argument
hrc_EUR_qc/hrc_EUR_qc.prune_rand2M_rep[1-20].snps           # ``--extract`` argument
```

Additionally, [here](https://github.com/precimed/mixer_private_docker/tree/main/reference/) one may download there are gene and gene-set definitions, derived as explained [here](https://github.com/ofrei/genesets/blob/main/prepare_genesets_v2.ipynb):

```
gsa-mixer-baseline-annot_10mar2023.csv           # ``--go-file`` (baseline model)
gsa-mixer-gene-annot_10mar2023.csv               # ``--go-file`` (model)
gsa-mixer-gene-annot_10mar2023-10kb-overlap.csv  # ``--go-file-overlap`` (genes overlap after 10kb expansion)
gsa-mixer-geneset-annot_10mar2023.csv            # ``--go-file-test`` (only gene-sets)
gsa-mixer-genesetLOO-annot_10mar2023.csv         # ``--go-file-test`  (only gene-sets, with leave-one-gene-out)
gsa-mixer-hybrid-annot_10mar2023.csv             # ``--go-file-test`  (genes and gene-sets)
gsa-mixer-hybridLOO-annot_10mar2023.csv          # ``--go-file-test`  (genes and gene-sets, with leave-one-gene-out)
magma-gene-annot_10mar2023.csv                   # gsa-mixer-gene-annot_10mar2023.csv converted to MAGMA format
magma-geneset-annot_10mar2023.csv                # gsa-mixer-geneset-annot_10mar2023.csv converted to MAGMA format
```

After downloading reference file we advice using ``--savelib-file`` option as shown below to generate ``.bin`` files,
with compressed representation of the reference. After that loading reference is possible with ``--loadlib-file``, providing considerable speedup over passing ``--ld-file`` argument.
In the following example the reference is saved to ``bin/`` folder - you need to change your own location.

The reference needs to be saved in two formats. The following example produces ``.bin`` files for ``plsa`` analysis,yielding its own ``.bin`` file for each chromosome. The ``--savelib-file`` argument must include ``@`` symbol which will be replaced with an actual chromosome label.
```
${MIXER_PY} plsa \
      --bim-file ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_chr@_EUR_qc.bim \
      --ld-file ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_chr@_EUR_qc.run1.ld \
      --use-complete-tag-indices --exclude-ranges [] \
      --savelib-file bin/1000G.EUR.QC.@.bin \
      --out bin/1000G.EUR.QC.@
```

The following example produces ``.bin`` file for ``fit1``,``fit2``,``test1``,``test2`` steps, yielding a single ``.bin`` file combined across all chromosomes. The ``--savelib-file`` argument does not need to include ``@`` symbol (if it does, the ``@`` symbol will stay unchanged, and simply be part of the output file name):
```
${MIXER_PY} fit1 \
      --bim-file ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_chr@_EUR_qc.bim \
      --ld-file ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_chr@_EUR_qc.run1.ld \
      --use-complete-tag-indices --exclude-ranges [] \
      --savelib-file bin/1000G.EUR.QC.@.bin \
      --out bin/1000G.EUR.QC.@
```

## (optional) Generate LD Reference

For real-data analysis reference data can usually be downloaded via the links provided above.
Otherwise GSA-MiXeR reference files can be prepared from plink bfile using ``mixer.py ld`` and ``mixer.py snps`` commands as described below. It's important that the reference genotypes contain unrelated individuals only, constrained to a single population.
Note that ``@`` symbol must remain as it is in all commands, i.e. you don't need to exchange it with a specific chromosome label.

Compute LD matrices (one per chromosome), later to be used with ``--ld-file`` argument.
```
#!/bin/bash
#SBATCH --job-name=gsamixer
#SBATCH --account=p697
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --cpus-per-task=8
#SBATCH --array=1-22

export CHR=${SLURM_ARRAY_TASK_ID}
export MIXER_SIF=mixer.sif
export MIXER_PY="singularity exec --home $PWD:/home ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"

${MIXER_PY} ld --bfile chr${CHR} --r2min 0.01 --ldscore-r2min 0.0001 --ld-window-kb 10000 --out chr${CHR}.ld
```

Compute SNPs subsets (one for each of 20 random repeats), later to be used with ``--extract`` argument.
```
#!/bin/bash
#SBATCH --job-name=gsamixer
#SBATCH --account=p697
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --cpus-per-task=8
#SBATCH --array=1-20

export REP=${SLURM_ARRAY_TASK_ID}
export MIXER_SIF=mixer.sif
export MIXER_PY="singularity exec --home $PWD:/home ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"

${MIXER_PY} snps --bim-file chr${CHR} --ld-file chr@ --chr2use 1-22 --maf 0.05 --subset 3000000 --seed $REP --out rep${REP}.snps
```

If only a random subset of SNPs is needed it's faster to use linux's ``cut`` and ``shuf`` commands:
```
for i in {1..20}
do 
  cat hrc_chr*_EUR_qc.bim | cut -f 2 | shuf | head -n 2000000 | sort > hrc_EUR_qc.prune_rand2M_rep$i.snps
done
```

Full command-line reference for ``mixer.py ld`` and ``mixer.py snps`` is as follows:
```
usage: mixer.py ld [-h] [--out OUT] [--lib LIB] [--log LOG] [--bfile BFILE]
                   [--r2min R2MIN] [--ldscore-r2min LDSCORE_R2MIN]
                   [--ld-window-kb LD_WINDOW_KB] [--ld-window LD_WINDOW]

optional arguments:
  -h, --help            show this help message and exit
  --out OUT             prefix for the output files, such as <out>.json
                        (default: mixer);
  --lib LIB             path to libbgmg.so plugin (default: libbgmg.so); can
                        be also specified via BGMG_SHARED_LIBRARY env
                        variable.
  --log LOG             file to output log (default: <out>.log); NB! if --log
                        points to an existing file the new lines will be
                        appended to it at the end of the file.
  --bfile BFILE         path to plink bfile (required argument)
  --r2min R2MIN         r2 values above this threshold will be stored in
                        sparse LD format (default: 0.05)
  --ldscore-r2min LDSCORE_R2MIN
                        r2 values above this threshold (and below --r2min)
                        will be stored as LD scores that contribute to the
                        cost function via an infinitesimal model (default:
                        0.001)
  --ld-window-kb LD_WINDOW_KB
                        limit window similar to --ld-window-kb in 'plink r2';
                        0 will disable this constraint (default: 0); either
                        ld-window-kb or ld-window argument must be provided
  --ld-window LD_WINDOW
                        limit window similar to --ld-window in 'plink r2'; 0
                        will disable this constraint (default: 0); either ld-
                        window-kb or ld-window argument must be provided
```

```
usage: mixer.py snps [-h] [--out OUT] [--lib LIB] [--log LOG]
                     [--bim-file BIM_FILE] [--ld-file LD_FILE]
                     [--chr2use CHR2USE] [--r2 R2] [--maf MAF]
                     [--subset SUBSET] [--seed SEED]

optional arguments:
  -h, --help           show this help message and exit
  --out OUT            prefix for the output files, such as <out>.json
                       (default: mixer);
  --lib LIB            path to libbgmg.so plugin (default: libbgmg.so); can be
                       also specified via BGMG_SHARED_LIBRARY env variable.
  --log LOG            file to output log (default: <out>.log); NB! if --log
                       points to an existing file the new lines will be
                       appended to it at the end of the file.
  --bim-file BIM_FILE  plink bim file (required argument); defines the
                       reference set of SNPs used for the analysis. Marker
                       names must not have duplicated entries. May contain
                       symbol '@', which will be replaced by an actual
                       chromosome label.
  --ld-file LD_FILE    file with linkage disequilibrium information, generated
                       via 'mixer.py ld' command (required argument); may
                       contain symbol '@', similarly to --bim-file argument.
  --chr2use CHR2USE    chromosome labels to use (default: 1-22); chromosome
                       must be labeled by integer, i.e. X and Y are not
                       acceptable; example of valid arguments: '1,2,3' or
                       '1-4,12,16-20'
  --r2 R2              LD r2 threshold for prunning SNPs (default: 0.8)
  --maf MAF            minor allele frequence (MAF) threshold (default: 0.05)
  --subset SUBSET      number of SNPs to randomly select (default: 2000000)
  --seed SEED          random seed (default: None)
```

## Perform GSA-MiXeR and MAGMA analyses

First step is to split summary statistics into one file per chromosome, as follows:

```
export MIXER_SIF=/ess/p697/data/durable/s3-api/github/norment/ofrei_repo/2023_03_27/mixer.sif
export MIXER_PY="singularity exec --home $PWD:/home --bind /ess/p697:/ess/p697 ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"
export SUMSTATS_FOLDER=/ess/p697/cluster/users/ofrei/ukbdata/projects/plsa_mixer/sumstats
export SUMSTATS_FILE=PGC_SCZ_0518_EUR
#export SUMSTATS_FILE=GIANT_HEIGHT_2018_UKB

${MIXER_PY} split_sumstats \
    --trait1-file ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.sumstats.gz
    --out ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.chr@.sumstats.gz
```

Then following two SLURM scripts will apply GSA-MiXeR and MAGMA analysis.
Commands below use reference data described in [Download LD Reference Files](#download-ld-reference-files) section.
MAGMA software is also included in ``mixer.sif`` container.

[scripts/GSA_MIXER.job](scripts/GSA_MIXER.job):
```
#!/bin/bash
#SBATCH --job-name=plsasimu
#SBATCH --account=p697_norment
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000M  # 502 GB available => 502*1024/64 = 8032 MB max per core
#SBATCH --cpus-per-task=8    # remember to update --threads argument below!
#SBATCH --array=1-20
##SBATCH --array=1

source /cluster/bin/jobsetup

test $SCRATCH && module load singularity/3.7.1 

export SINGULARITY_BIND=
export MIXER_SIF=/ess/p697/data/durable/s3-api/github/norment/ofrei_repo/2023_03_27/mixer.sif
md5sum ${MIXER_SIF}  # b99f7b027b0448d9338f9506bac09c66

export MIXER_PY="singularity exec --home $PWD:/home --bind /ess/p697:/ess/p697 ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"

export SUMSTATS_FOLDER=/ess/p697/cluster/users/ofrei/ukbdata/projects/plsa_mixer/sumstats

export SUMSTATS_FILE=PGC_SCZ_0518_EUR
#export SUMSTATS_FILE=GIANT_HEIGHT_2018_UKB

export OUT_FOLDER=/ess/p697/cluster/users/ofrei/ukbdata/projects/plsa_mixer/real_jan22_main_run1
export REFERENCE_FOLDER=/ess/p697/data/durable/s3-api/github/precimed/mixer_private_docker/reference
export REP=${SLURM_ARRAY_TASK_ID}

# just to test if commands run, add the following:
# --chr2use 21-22 --adam-epoch 3 3 --adam-step 0.064 0.032

# baseline
${MIXER_PY} plsa \
        --trait1-file ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.chr@.sumstats.gz \
        --out ${OUT_FOLDER}/${SUMSTATS_FILE}_baseline_rep${REP} \
        --bim-file ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_chr@_EUR_qc.bim \
        --ld-file ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_chr@_EUR_qc.run1.ld \
        --extract ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_EUR_qc.prune_rand2M_rep${REP}.snps \
        --go-file ${REFERENCE_FOLDER}/gsa-mixer-baseline-annot_27jan2022.csv \
        --go-file-test ${REFERENCE_FOLDER}/gsa-mixer-hybridLOO-annot_27jan2022.csv \
        --annot-file ${REFERENCE_FOLDER}/hrc_EUR_qc/baseline_v2.2_hrc_chr@_EUR_qc.annot.gz \
        --annot-file-test ${REFERENCE_FOLDER}/hrc_EUR_qc/baseline_v2.2_hrc_chr@_EUR_qc.annot.gz \
        --z1max 9.336 --sig2-zeroL 0 --s-value -0.25 --pi-value 1.0 --l-init -0.25 \
        --seed $((REP+1000)) --threads 8

# enrichment model
${MIXER_PY} plsa \
        --trait1-file ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.chr@.sumstats.gz \
        --out ${OUT_FOLDER}/${SUMSTATS_FILE}_model_rep${REP} \
        --bim-file ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_chr@_EUR_qc.bim \
        --ld-file ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_chr@_EUR_qc.run1.ld \
        --extract ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_EUR_qc.prune_rand2M_rep${REP}.snps \
        --go-file ${REFERENCE_FOLDER}/gsa-mixer-gene-annot_27jan2022.csv \
        --go-file-test ${REFERENCE_FOLDER}/gsa-mixer-hybridLOO-annot_27jan2022.csv \
        --annot-file-test ${REFERENCE_FOLDER}/hrc_EUR_qc/baseline_v2.2_hrc_chr@_EUR_qc.annot.gz \
        --z1max 9.336 --sig2-zeroL 0 --s-value -0.25 --pi-value 1.0 --l-init -0.25 \
        --seed $((REP+1000)) --threads 8
```

[scripts/MAGMA.job](scripts/MAGMA.job):
```
#!/bin/bash
#SBATCH --job-name=magma
#SBATCH --account=p697_norment
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000M  # 502 GB available => 502*1024/64 = 8032 MB max per core
#SBATCH --cpus-per-task=4

source /cluster/bin/jobsetup

test $SCRATCH && module load singularity/3.7.1 

export SINGULARITY_BIND=
export MIXER_SIF=/ess/p697/data/durable/s3-api/github/norment/ofrei_repo/2023_03_27/mixer.sif
md5sum ${MIXER_SIF}  # b99f7b027b0448d9338f9506bac09c66

export MAGMA="singularity exec --home $PWD:/home --bind /ess/p697:/ess/p697 ${MIXER_SIF} magma"

export SUMSTATS_FOLDER=/ess/p697/cluster/users/ofrei/ukbdata/projects/plsa_mixer/sumstats

export SUMSTATS_FILE=PGC_SCZ_0518_EUR
#export SUMSTATS_FILE=GIANT_HEIGHT_2018_UKB

export OUT_FOLDER=/ess/p697/cluster/users/ofrei/ukbdata/projects/plsa_mixer/real_jan22_magma_run1

export REFERENCE_FOLDER=/ess/p697/data/durable/s3-api/github/precimed/mixer_private_docker/reference
export MAGMA_BFILE=/ess/p697/cluster/projects/moba_qc_imputation/resources/HRC/plink/hrc_EUR_qc_keep1k
export MAGMA_GENE_LOC=${REFERENCE_FOLDER}/magma-gene-annot_27jan2022.csv
export MAGMA_SET_ANNOT=${REFERENCE_FOLDER}/magma-geneset-annot_27jan2022.csv
export OUT=${OUT_FOLDER}/${SUMSTATS_FILE}

zcat  ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.sumstats.gz >  ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.sumstats

$MAGMA --snp-loc ${MAGMA_BFILE}.bim \
       --gene-loc ${MAGMA_GENE_LOC} \
       --out $OUT.magma.step1 \
       --annotate window=10 
$MAGMA --pval ${SUMSTATS_FOLDER}/${SUMSTATS_FILE}.sumstats snp-id=SNP pval=PVAL ncol=N \
       --bfile ${MAGMA_BFILE} \
       --gene-annot $OUT.magma.step1.genes.annot \
       --out $OUT.magma.step2
$MAGMA --gene-results $OUT.magma.step2.genes.raw \
       --set-annot ${MAGMA_SET_ANNOT} \
       --out $OUT.magma
```

## Output formats of GSA-MiXeR analysis

The above scripts produce the following output files:

```
PGC_SCZ_0518_EUR_baseline_rep[1-20].go_test_enrich.csv
PGC_SCZ_0518_EUR_model_rep[1-20].go_test_enrich.csv
```

Final step of the analysis is to compute ratios between ``model`` and ``baseline``
using [scripts/gsa_mixer_combine.py](scripts/gsa_mixer_combine.py) (not integrated into .sif container).
Later this will be available as follows:
```
export SUMSTATS_FILE=PGC_SCZ_0518_EUR
export MIXER_FIGURES_PY="singularity exec --home $PWD:/home ${MIXER_SIF} python /tools/mixer/precimed/mixer_figures.py"
MIXER_FIGURES_PY combine --go-file-out-baseline ${SUMSTATS_FILE}_baseline_rep@.go_test_enrich.csv \
                         --go-file-out-model ${SUMSTATS_FILE}_model_rep@.go_test_enrich.csv \
                         --out ${SUMSTATS_FILE}.go_test_enrich.csv
```

## Command-line reference

```
usage: mixer.py plsa [-h] [--out OUT] [--lib LIB] [--log LOG]
                     [--bim-file BIM_FILE] [--ld-file LD_FILE]
                     [--chr2use CHR2USE] [--extract EXTRACT]
                     [--exclude EXCLUDE]
                     [--exclude-ranges EXCLUDE_RANGES [EXCLUDE_RANGES ...]]
                     [--allow-ambiguous-snps] [--trait1-file TRAIT1_FILE]
                     [--z1max Z1MAX] [--randprune-n RANDPRUNE_N]
                     [--randprune-r2 RANDPRUNE_R2]
                     [--randprune-maf RANDPRUNE_MAF] [--seed SEED]
                     [--threads THREADS [THREADS ...]]
                     [--kmax KMAX [KMAX ...]] [--annot-file ANNOT_FILE]
                     [--annot-file-test ANNOT_FILE_TEST] [--go-file GO_FILE]
                     [--go-file-test GO_FILE_TEST]
                     [--go-all-genes-label GO_ALL_GENES_LABEL]
                     [--go-extend-bp GO_EXTEND_BP] [--sig2-zeroL SIG2_ZEROL]
                     [--s-value S_VALUE] [--l-value L_VALUE]
                     [--pi-value PI_VALUE] [--s-init S_INIT] [--l-init L_INIT]
                     [--pi-init PI_INIT] [--annot-p ANNOT_P] [--gene-p GENE_P]
                     [--adam-epoch ADAM_EPOCH [ADAM_EPOCH ...]]
                     [--adam-beta1 ADAM_BETA1] [--adam-beta2 ADAM_BETA2]
                     [--adam-eps ADAM_EPS]
                     [--adam-step ADAM_STEP [ADAM_STEP ...]] [--adam-disable]
                     [--load-params-file LOAD_PARAMS_FILE] [--make-snps-file]

optional arguments:
  -h, --help            show this help message and exit
  --out OUT             prefix for the output files, such as <out>.json
                        (default: mixer);
  --lib LIB             path to libbgmg.so plugin (default: libbgmg.so); can
                        be also specified via BGMG_SHARED_LIBRARY env
                        variable.
  --log LOG             file to output log (default: <out>.log); NB! if --log
                        points to an existing file the new lines will be
                        appended to it at the end of the file.
  --bim-file BIM_FILE   plink bim file (required argument); defines the
                        reference set of SNPs used for the analysis. Marker
                        names must not have duplicated entries. May contain
                        symbol '@', which will be replaced by an actual
                        chromosome label.
  --ld-file LD_FILE     file with linkage disequilibrium information,
                        generated via 'mixer.py ld' command (required
                        argument); may contain symbol '@', similarly to --bim-
                        file argument.
  --chr2use CHR2USE     chromosome labels to use (default: 1-22); chromosome
                        must be labeled by integer, i.e. X and Y are not
                        acceptable; example of valid arguments: '1,2,3' or
                        '1-4,12,16-20'
  --extract EXTRACT     (optional) File with variants to include in the
                        analysis. By default, all variants are included. This
                        applies to GWAS tag SNPs, however LD is still computed
                        towards the entire reference provided by --bim-file.
                        This applies before --exclude, so a variant listed in
                        --exclude won't be re-introduced if it's also present
                        in --extract list.
  --exclude EXCLUDE     (optional) File with variants to exclude from the
                        analysis.
  --exclude-ranges EXCLUDE_RANGES [EXCLUDE_RANGES ...]
                        (default: ['MHC']) exclude SNPs in ranges of base pair
                        position; the syntax is chr:from-to, for example
                        6:25000000-35000000; multiple regions can be excluded;
                        "chr" prefix prior to chromosome label, as well as KB
                        and MB suffices are allowed, e.g. chr6:25-35MB is a
                        valid exclusion range. Some special case regions are
                        also supported, for example "--exclude-ranges MHC
                        APOE".To overwrite the defaul, pass "--exclude ranges
                        []".
  --allow-ambiguous-snps
                        advanced option (expert use only); a flag allowing to
                        include A/T and C/G SNPs in fit procedure.
  --trait1-file TRAIT1_FILE
                        GWAS summary statistics for the first trait (required
                        argument); for 'plsa' analysis it is recommended to
                        split GWAS summary statistics per chromosome; if this
                        is done then --trait1-file should contain symbol '@',
                        which will be replaced by an actual chromosome label.
  --z1max Z1MAX         right-censoring threshold for the first trait
                        (default: None); recommended setting: '--z1max 9.336'
                        (equivalent to p-value 1e-20)
  --randprune-n RANDPRUNE_N
                        number of random pruning iterations (default: 64)
  --randprune-r2 RANDPRUNE_R2
                        threshold for random pruning (default: 0.1)
  --randprune-maf RANDPRUNE_MAF
                        threshold for minor allele frequency (default: 0.05);
                        applies to tag SNPs to include in the fit procedure
  --seed SEED           random seed (default: None).
  --threads THREADS [THREADS ...]
                        specify how many concurrent threads to use for
                        computations; (default: total number of CPU cores
                        available)
  --kmax KMAX [KMAX ...]
                        number of sampling iterations for log-likelihod and
                        posterior delta (default: 20000)
  --annot-file ANNOT_FILE
                        (optional) path to binary annotations in LD score
                        regression format, i.e. <path>/baseline.@.annot.gz for
                        fitting enrichment model model. This must include the
                        first column with all ones ('base' annotation category
                        covering the entire genome). Enrichment scores
                        computed for --annot-file will be saved to
                        <out>.annot_enrich.csv
  --annot-file-test ANNOT_FILE_TEST
                        (optional) path to binary annotations in LD score
                        regression format, i.e. <path>/baseline.@.annot.gz for
                        evaluating enrichment. If provided, enrichment scores
                        computed for --annot-file-test will be saved to
                        <out>.annot_test_enrich.csv.
  --go-file GO_FILE     (optional) path to GO antology file for fitting
                        enrichment model model. The format is described in the
                        documentation. 'base' category that covers entire
                        genome will be added automatically. Enrichment scores
                        computed for --go-file will be saved to
                        <out>.go_enrich.csv
  --go-file-test GO_FILE_TEST
                        (optional) path to an additional GO antology file for
                        evaluating enrichment, same convention as --go-file.
                        If provided, enrichment scores computed for --go-test-
                        file will be saved to <out>.go_test_enrich.csv
  --go-all-genes-label GO_ALL_GENES_LABEL
                        reference gene-set to calibrate fold enrichment, e.g.
                        allowing to compute enrichment w.r.t. the set of all
                        coding genes (default:'base')
  --go-extend-bp GO_EXTEND_BP
                        extends each gene by this many base pairs, defining a
                        symmetric window up and downstream (default: 10000)
  --sig2-zeroL SIG2_ZEROL
                        (optional) constraint for 'sig2_zeroL' parameter;
                        recommended setting for gene-set enrichment analysis:
                        '--sig2-zeroL 0'
  --s-value S_VALUE     (optional) constraint for the 's' parameter;
                        recommended setting for gene-set enrichment analysis:
                        '--s-value -0.25'
  --l-value L_VALUE     (optional) constraint for the 'l' parameter
  --pi-value PI_VALUE   (optional) constraint for the 'pi' parameter;
                        recommended setting for gene-set enrichment analysis:
                        '--pi-value 1.0'
  --s-init S_INIT       initial value for the 's' parameter (default: -0.25);
                        does not apply if --s-value is provided)
  --l-init L_INIT       initial value for the 'l' parameter (default: -0.125);
                        does not apply if --l-value is provided)
  --pi-init PI_INIT     initial value for the 'pi' parameter (default: 0.001);
                        does not apply if --pi-value is provided)
  --annot-p ANNOT_P     power factor for sigma2 aggregation in overlapping
                        annotations (default: 1)
  --gene-p GENE_P       power factor for sigma2 aggregation in overlapping
                        gene-sets (default: 1)
  --adam-epoch ADAM_EPOCH [ADAM_EPOCH ...]
                        number of iterations in ADAM procedure (default: [10,
                        10, 10, 10, 10, 10, 10, 10, 10, 10])
  --adam-beta1 ADAM_BETA1
                        beta_1 parameter in ADAM procedure (default: 0.9)
  --adam-beta2 ADAM_BETA2
                        beta_2 parameter in ADAM procedure (default: 0.99)
  --adam-eps ADAM_EPS   epsilon parameter in ADAM procedure (default: 1e-08)
  --adam-step ADAM_STEP [ADAM_STEP ...]
                        step parameter in ADAM procedure (default: [0.064,
                        0.032, 0.016, 0.008, 0.004, 0.002, 0.001, 0.0005,
                        0.00025, 0.0001])
  --adam-disable        a flag allowing to disable optimization; typical
                        usecase would be in conjunction with these flags: '--
                        adam-disable --load-params-file <out-of-a-previous-
                        run>.json --make-snps-file --allow-ambiguous-snps'
  --load-params-file LOAD_PARAMS_FILE
                        (optional) params of the fitted model to load;
                        expected to be from a 'mixer.py plsa' run
  --make-snps-file      a flag allowing to generate file with per-SNP
                        estimates; will generate <out>.snps.csv output file
```
