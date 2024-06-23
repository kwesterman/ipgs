#!/bin/bash

#$ -l h_vmem=80G
#$ -l h_rt=6:00:00
#$ -o ../reports

#$ -j y
#$ -cwd



pheno=$1


MAF=0.005
phenoFile=../data/processed/ukb_training_set.csv


vgwas_dir=../data/processed/vgwas
scratch=/broad/hptmp/gervis

fn_prefix=${vgwas_dir}/${pheno}



source /broad/software/scripts/useuse
use R-4.1

mkdir -p ../reports




# prep pheno & covar file for QUAIL
Rscript --vanilla vgwas/prep_vgwas.R ${pheno} ${phenoFile}


# Calculate integrated rank score using QUAIL
Rscript ./QUAIL/Step1_QUAIL_rank_score.R \
--pheno ${fn_prefix}/_pheno.txt \
--covar ${fn_prefix}_covars.txt \
--output ${fn_prefix}_rankscore.txt \
--num_levels 2000 \
--num_cores 5




#EOF
