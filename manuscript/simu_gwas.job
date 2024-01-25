#!/bin/bash
#SBATCH --job-name=simugwa5
#SBATCH --account=p697
##SBATCH --account=p697
#SBATCH --time=12:00:00  
##SBATCH --time=4:00:00  
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000M  # 502 GB available => 502*1024/64 = 8032 MB max per core
#SBATCH --cpus-per-task=16
##SBATCH --array 1-200
##SBATCH --array 1-150
##SBATCH --array 1-400
##SBATCH --array 1-1000
#SBATCH --array 1-200
##SBATCH --array 1-200
##SBATCH --array 1-50

export THREADS=16

source /cluster/bin/jobsetup
test $SCRATCH && module load singularity/3.7.1 

export FOLDER="/ess/p697/cluster/users/ofrei/2023_02_06_GSA_MiXeR_natgen_revisions/simu_snps"

# 200 runs to run10 folder; HSQ=0.1, 0.4, 0.7
export TASKLIST="/ess/p697/cluster/users/ofrei/gsa_mixer_analysis/task_list_run10a.txt"

# 200 runs in run10b folders; HSQ=0.1, 0.4, 0.7
#export TASKLIST="/ess/p697/cluster/users/ofrei/gsa_mixer_analysis/task_list_run10b.txt"

# 150 runs to run11 folder; HSQ=0.7
#export TASKLIST="/ess/p697/cluster/users/ofrei/gsa_mixer_analysis/task_list_run11a.txt"

# 400 runs to run13 folder; HSQ=0.7
#export TASKLIST="/ess/p697/cluster/users/ofrei/gsa_mixer_analysis/task_list_run13.txt"

# 1000 runs to run14 folder; HSQ=0.7
#export TASKLIST="/ess/p697/cluster/users/ofrei/gsa_mixer_analysis/task_list_run14a.txt"
#export TASKLIST="/ess/p697/cluster/users/ofrei/gsa_mixer_analysis/task_list_run14b.txt"
#export TASKLIST="/ess/p697/cluster/users/ofrei/gsa_mixer_analysis/task_list_run14c.txt"

export FNAME=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${TASKLIST})
echo "processing line ${SLURM_ARRAY_TASK_ID}: ${FNAME}... "

#export HSQ="0.1"
#export HSQ="0.4"
#export HSQ="0.7"

export HSQ=$1
echo "HSQ=${HSQ}"

export TMP_GWAS_OUT="${SCRATCH}/${FNAME}_hsq=${HSQ}"
export FNAME_GWAS_PREFIX="${FOLDER}/gwas_run14/${FNAME}"
export FNAME_GSA_BASE_OUT="${FOLDER}/base_run14/${FNAME}_hsq=${HSQ}"
export FNAME_GSA_FULL_OUT="${FOLDER}/full_run14/${FNAME}_hsq=${HSQ}"
export FNAME_GSA_HESS_OUT="${FOLDER}/hess_run14/${FNAME}_hsq=${HSQ}"

export FNAME_MAGMA_OUT="${FOLDER}/magma_run11_null/${FNAME}_hsq=${HSQ}"

export FNAME_GWAS_OUT="${FNAME_GWAS_PREFIX}_hsq=${HSQ}" 

export SUMSTATS_FILE=${FNAME_GWAS_OUT}.chr@.sumstats.gz
export SUMSTATS_MAGMA_FILE=${FNAME_GWAS_OUT}.chr@.sumstats
export SEED=123

export CAUSAL_VARIANTS_FOLDER=${FOLDER}/causal_variants3
export CAUSALS_BETA="${CAUSAL_VARIANTS_FOLDER}/${FNAME}.causals_beta.csv"
export GO_TEST_FILE="${CAUSAL_VARIANTS_FOLDER}/${FNAME}.genes.csv"
export MAGMA_SET_ANNOT="${CAUSAL_VARIANTS_FOLDER}/${FNAME}.magma.csv"

export MIXER_SIF="/ess/p697/data/durable/s3-api/github/norment/ofrei_repo/2023_03_27/mixer.sif"
export SIMU_LINUX="singularity exec --home $PWD:/home ${MIXER_SIF} simu_linux"
export PLINK2="singularity exec --home $PWD:/home ${MIXER_SIF} plink2"
export PYTHON="singularity exec --home $PWD:/home ${MIXER_SIF} python"
export MAGMA="singularity exec --home $PWD:/home --bind /ess/p697:/ess/p697 ${MIXER_SIF} magma"
#export MIXER_PY="singularity exec --home $PWD:/home  ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"
#export MIXER_PY="singularity exec --home $PWD:/home ${MIXER_SIF} python /ess/p697/data/durable/s3-api/github/precimed/mixer_private/precimed/mixer.py"
export MIXER_PY="singularity exec --home $PWD:/home ${MIXER_SIF} python /ess/p697/cluster/users/ofrei/2023_02_06_GSA_MiXeR_natgen_revisions/precimed/mixer.py"

# perform GWAS
if [ ! -f "${FNAME_GWAS_PREFIX}_hsq=0.1.pheno" ]; then
export BFILE="/ess/p697/cluster/users/ofrei/ukbdata/projects/plsa_mixer/ukb_genetics_qc/ukb_bed/ukb_imp_chr@_v3_qc"
${SIMU_LINUX} --seed $SEED --bfile-chr $BFILE --qt --causal-variants ${CAUSALS_BETA} --hsq ${HSQ} --out ${FNAME_GWAS_OUT}
#${PYTHON} computed_pheno_from_genepheno_v3.py ${CAUSAL_VARIANTS_FOLDER}/${FNAME} ${FNAME_GWAS_PREFIX}
else
  echo "PHENO output already exists, skip (${FNAME_GWAS_PREFIX}_hsq=XXX.pheno)"
fi

if [ ! -f "${SUMSTATS_MAGMA_FILE}" ]; then
export PLINK_EXTRACT="/cluster/projects/p697/users/ofrei/2023_02_06_GSA_MiXeR_natgen_revisions/simu_snps/base_maf=0.05_LDr2=0.90.snps"
{ awk '{print $1}' ${CAUSALS_BETA}; cat ${PLINK_EXTRACT}; } | sort | uniq > ${TMP_GWAS_OUT}.gwas.snps
for chri in {1..22};
do
  export BFILE="/ess/p697/cluster/users/ofrei/ukbdata/projects/plsa_mixer/ukb_genetics_qc/ukb_bed/ukb_imp_chr${chri}_v3_qc"
  ${PLINK2} --glm allow-no-covars --extract ${TMP_GWAS_OUT}.gwas.snps --bfile $BFILE --pheno ${FNAME_GWAS_OUT}.pheno  --pheno-name trait1 --out ${TMP_GWAS_OUT}.chr${chri}
done 

$PYTHON concat_plink_and_convert.py ${TMP_GWAS_OUT}.chr@.trait1.glm.linear ${FNAME_GWAS_OUT}.chr@

for chri in {1..22};
do
  echo "purge temp files for chr${chri}"
  rm ${TMP_GWAS_OUT}.chr${chri}.trait1.glm.linear
  rm ${TMP_GWAS_OUT}.chr${chri}.log
done 
rm ${TMP_GWAS_OUT}.gwas.snps

else
  echo "GWAS output already exists, skip (${SUMSTATS_MAGMA_FILE})"
fi  # perform GWAS

export REFERENCE_FOLDER=/ess/p697/data/durable/s3-api/github/precimed/mixer_private_docker/reference
export ANNOT_VERSION=annot_10mar2023
export EXTRA_FLAGS=" --seed 1001  --exclude-ranges MHC"

# just to test if commands run, add the following:
#export EXTRA_JUNK="--chr2use 20-22 --adam-epoch 3 3 --adam-step 0.064 0.032"
export EXTRA_JUNK=

# hrc-based panel
#export COMMON_ARGS=""
#export COMMON_ARGS="${COMMON_ARGS} --bim-file ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_chr@_EUR_qc.bim "
#export COMMON_ARGS="${COMMON_ARGS} --ld-file ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_chr@_EUR_qc.run1.ld "
#export COMMON_ARGS="${COMMON_ARGS} --use-complete-tag-indices --loadlib-file /ess/p697/cluster/projects/moba_qc_imputation/resources/HRC/plink/tsd_libfile_hrc_chr@.bin "
#export COMMON_ARGS="${COMMON_ARGS} --annot-file ${REFERENCE_FOLDER}/hrc_EUR_qc/baseline_v2.2_hrc_chr@_EUR_qc.annot.gz "
#export COMMON_ARGS="${COMMON_ARGS} --annot-file-test ${REFERENCE_FOLDER}/hrc_EUR_qc/baseline_v2.2_hrc_chr@_EUR_qc.annot.gz "

# ukb-based panel
export COMMON_ARGS=""
export COMMON_ARGS="${COMMON_ARGS} --bim-file ${REFERENCE_FOLDER}/ukb_EUR_qc/ukb_imp_chr@_v3_qc.bim "
export COMMON_ARGS="${COMMON_ARGS} --ld-file ${REFERENCE_FOLDER}/ukb_EUR_qc/ukb_imp_chr@_v3_qc.run1.ld "
export COMMON_ARGS="${COMMON_ARGS} --use-complete-tag-indices --loadlib-file /ess/p697/cluster/users/ofrei/ukbdata/projects/plsa_mixer/lib_bin_randprune64/plsa_ukb_libbgmg_state_full_chr@.bin "

#export ANNOT_FILE="${REFERENCE_FOLDER}/ukb_EUR_qc/baseline_v2.2_ukb_imp_chr@_v3_qc.annot.gz"
export ANNOT_FILE="${REFERENCE_FOLDER}/ukb_EUR_qc/MAFbins_ukb_imp_chr@_v3_qc.annot.gz"
#export ANNOT_FILE="${REFERENCE_FOLDER}/ukb_EUR_qc/MAFbins_LLDbins_ukb_imp_chr@_v3_qc.annot.gz"

export COMMON_ARGS="${COMMON_ARGS} --annot-file ${ANNOT_FILE} --annot-file-test ${ANNOT_FILE} "

# fit GSA BASE MiXeR model

if [ ! -f "${FNAME_GSA_BASE_OUT}_baseline.json" ]; then

#${MIXER_PY} plsa --gsa-base wth --l-value 0 --s-value 0 \

${MIXER_PY} plsa --fit sig2-zeroA annot gene --l-value 0 --s-value 0  --h2-init-calibration 0.1 \
        --trait1-file ${SUMSTATS_FILE} \
        --out ${FNAME_GSA_BASE_OUT}_baseline \
        --go-file ${REFERENCE_FOLDER}/gsa-mixer-baseline-${ANNOT_VERSION}.csv \
        --go-file-test ${GO_TEST_FILE} \
        --threads ${THREADS} ${EXTRA_FLAGS} ${EXTRA_JUNK} ${COMMON_ARGS} \
        --se-samples 0

else
  echo "GSA BASE output already exists, skip (${FNAME_GSA_BASE_OUT}_baseline.json)"
fi

# fit GSA FULL MiXeR model
if [ ! -f "${FNAME_GSA_FULL_OUT}_model.json" ]; then

# --gsa-full excluding --calc-loglike-diff-go (i.e. speed up due to not computing AICs)

${MIXER_PY} plsa --fit gene --constrain-base-gene-category --nullify-go-all-genes-sig2-gene --adam-separate-by-chromosome --adam-beta1 0.1 --adam-beta2 0.8 \
        --trait1-file ${SUMSTATS_FILE} \
        --out ${FNAME_GSA_FULL_OUT}_model \
        --go-file ${REFERENCE_FOLDER}/gsa-mixer-gene-${ANNOT_VERSION}.csv \
        --go-file-test ${GO_TEST_FILE} \
        --load-params-file ${FNAME_GSA_BASE_OUT}_baseline.json \
        --threads ${THREADS} ${EXTRA_FLAGS} ${EXTRA_JUNK} ${COMMON_ARGS} \
  
else
  echo "GSA FULL output already exists, skip (${FNAME_GSA_FULL_OUT}_model.json)"
fi  # fit GSA MiXeR model


# compute hessian and error bars
if [ ! -f "${FNAME_GSA_HESS_OUT}_model.go_test_enrich.csv" ]; then

${MIXER_PY} plsa --fit gene --constrain-base-gene-category --nullify-go-all-genes-sig2-gene --adam-separate-by-chromosome --adam-beta1 0.1 --adam-beta2 0.8 \
        --trait1-file ${SUMSTATS_FILE} \
        --out ${FNAME_GSA_HESS_OUT}_model \
        --go-file ${REFERENCE_FOLDER}/gsa-mixer-gene-${ANNOT_VERSION}.csv \
        --go-file-test ${GO_TEST_FILE} \
        --load-params-file ${FNAME_GSA_FULL_OUT}_model.json \
        --load-baseline-params-file ${FNAME_GSA_BASE_OUT}_baseline.json \
        --threads ${THREADS} ${EXTRA_FLAGS} ${EXTRA_JUNK} ${COMMON_ARGS} \
        --adam-disable --weights-file none --hardprune-r2 0.6 \

else
  echo "GSA HESS output already exists, skip (${FNAME_GSA_HESS_OUT}_model.go_test_enrich.csv)"
fi  # fit GSA MiXeR model


#export MAGMA_BFILE=/ess/p697/cluster/projects/moba_qc_imputation/resources/HRC/plink/hrc_EUR_qc_keep1k
export MAGMA_BFILE=/ess/p697/cluster/users/ofrei/ukbdata/projects/plsa_mixer/ukb_genetics_qc/ukb_bed/ukb_imp_chr@_v3_qc_keep1k
export MAGMA_GENE_LOC=${REFERENCE_FOLDER}/magma-gene-${ANNOT_VERSION}.csv

# MAGMA 
if [ ! -f "${FNAME_MAGMA_OUT}.magma.gsa.out" ]; then

$MAGMA --snp-loc ${MAGMA_BFILE}.bim \
       --gene-loc ${MAGMA_GENE_LOC} \
       --out ${FNAME_MAGMA_OUT}.magma.step1 \
       --annotate window=10 
$MAGMA --pval ${SUMSTATS_MAGMA_FILE} snp-id=SNP pval=PVAL ncol=N \
       --bfile ${MAGMA_BFILE} \
       --gene-annot ${FNAME_MAGMA_OUT}.magma.step1.genes.annot \
       --out ${FNAME_MAGMA_OUT}.magma.step2
$MAGMA --gene-results ${FNAME_MAGMA_OUT}.magma.step2.genes.raw \
       --set-annot ${MAGMA_SET_ANNOT} \
       --out ${FNAME_MAGMA_OUT}.magma
else
       echo "MAGMA output already exists, skip (${FNAME_MAGMA_OUT}.magma.gsa.out)"
fi # MAGMA

