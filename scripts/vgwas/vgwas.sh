#!/bin/bash

#$ -l h_vmem=80G
#$ -l h_rt=8:00:00
#$ -o ../reports


#$ -j y
#$ -cwd



CHR=$SGE_TASK_ID

MAF=$1
pheno=$2


opt=../opt
scratch=/broad/hptmp/gervis

source /broad/software/scripts/useuse
use R-4.1



# Run QUAIL Step 2 – vGWAS 
Rscript vgwas/QUAIL/Step2_QUAIL_vQTL_JEG.R \
--pheno ../data/vgwas/quail/ukb_${pheno}_quail.txt \
--pheno_rs ../data/vgwas/quail/ukb_${pheno}_rank_score.txt \
--covar ../data/vgwas/quail/ukb_covars_quail.txt \
--geno ${scratch}/plinkset/chr${CHR}_sel_maf${MAF} \
--freq ${scratch}/plinkset/chr${CHR}_sel_maf${MAF}.afreq \
--output ../data/vgwas/quail/ukb_chr${CHR}_${pheno}_rank_score \
--plink_path ${opt}/plink2


## EOF

