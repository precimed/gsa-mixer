#!/bin/bash
#SBATCH --job-name=gsamix
#SBATCH --account=p697_norment
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000M  # 502 GB available => 502*1024/64 = 8032 MB max per core
#SBATCH --cpus-per-task=8    # remember to update --threads argument below!
#SBATCH --array=1-36
##SBATCH --array=1-1

export THREADS=8

source /cluster/bin/jobsetup

test $SCRATCH && module load singularity/3.7.1 

export SUMSTATS_FOLDER=/ess/p697/cluster/users/ofrei/2023_02_06_GSA_MiXeR_natgen_revisions/sumstats_v3p1

export REFERENCE_FOLDER=/ess/p697/data/durable/s3-api/github/precimed/mixer_private_docker/reference
export ANNOT_VERSION=annot_10mar2023

export SINGULARITY_BIND=
export MIXER_SIF=/ess/p697/data/durable/s3-api/github/norment/ofrei_repo/2023_03_27/mixer.sif
export LDSC_SIF=/ess/p697/data/durable/s3-api/github/comorment/ldsc/containers/ldsc.sif

export MIXER_PY="singularity exec --home $PWD:/home --bind /ess/p697:/ess/p697 ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"
export PYTHON="singularity exec --home $PWD:/home --bind /ess/p697:/ess/p697 ${MIXER_SIF} python"

# /ess/p697/data/durable/s3-api/github/ofrei/ldsc/ldsc.py - custom version allowing --n-blocks 0 to disable jack-knife
export LDSC_PY="singularity exec --home $PWD:/home --bind /ess/p697:/ess/p697 ${LDSC_SIF} python /ess/p697/data/durable/s3-api/github/ofrei/ldsc/ldsc.py"
export MUNGE_PY="singularity exec --home $PWD:/home --bind /ess/p697:/ess/p697 ${LDSC_SIF} python /ess/p697/data/durable/s3-api/github/ofrei/ldsc/munge_sumstats.py"

#export JUPYTER="singularity exec --home $PWD:/home --bind /ess/p697:/ess/p697 ${MIXER_SIF} jupyter notebook --port 8989"
#export MIXER_PY="singularity exec --home $PWD:/home --bind /ess/p697:/ess/p697 ${MIXER_SIF} python /ess/p697/data/durable/s3-api/github/precimed/mixer_private/precimed/mixer.py"

export LOADLIB_FILE=/ess/p697/cluster/projects/moba_qc_imputation/resources/HRC/plink/tsd_libfile_hrc_chr@.bin

export ANNOT_FILE="${REFERENCE_FOLDER}/hrc_EUR_qc/baseline_v2.2_hrc_chr@_EUR_qc.annot.gz"


if false; then
${MIXER_PY} plsa --save-ldsc-reference \
        --bim-file ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_chr@_EUR_qc.bim \
        --ld-file ${REFERENCE_FOLDER}/hrc_EUR_qc/hrc_chr@_EUR_qc.run1.ld \
        --use-complete-tag-indices --loadlib-file ${LOADLIB_FILE} \
        --annot-file ${ANNOT_FILE} \
        --out baseline_v2.2_hrc_chr@_EUR_qc \
        --threads ${THREADS}
#        --go-file ldsc_scz_geneset${SLURM_ARRAY_TASK_ID}.csv \
fi

if false; then
${MUNGE_PY} --sumstats /ess/p697/cluster/users/ofrei/2023_02_06_GSA_MiXeR_natgen_revisions/sumstats_v3p1/PGC_SCZ_0518_EUR.sumstats \
--merge-alleles /ess/p697/data/durable/s3-api/github/comorment/ldsc/reference/w_hm3.snplist \
--out PGC_SCZ_0518_EUR --ignore B,OR --a1 EffectAllele --a2 OtherAllele --snp RSID
mv PGC_SCZ_0518_EUR.sumstats.gz PGC_SCZ_0518_EUR.sumstats
fi

export FRQFILE=/ess/p697/cluster/projects/moba_qc_imputation/resources/HRC/plink/hrc_chr@_EUR_qc
${LDSC_PY} --h2 PGC_SCZ_0518_EUR.sumstats \
	   --ref-ld-chr baseline_v2.2_hrc_chr@_EUR_qc,ldsc_scz_geneset${SLURM_ARRAY_TASK_ID}.chr@ \
	   --w-ld-chr /ess/p697/data/durable/s3-api/github/comorment/ldsc/reference/eur_w_ld_chr/ \
	   --overlap-annot --n-blocks 0 \
    	   --frqfile-chr ${FRQFILE} --out PGC_SCZ_0518_EUR_geneset${SLURM_ARRAY_TASK_ID}
