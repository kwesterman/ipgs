#!/bin/bash

#$ -l h_vmem=60G
#$ -l h_rt=24:00:00
#$ -o ../reports

#$ -j y
#$ -cwd



MAF=$1
phenoFile=$2
pheno=$3

opt=../opt

mkdir -p ../reports

source /broad/software/scripts/useuse
use R-4.1
use UGER


mkdir -p ../data/vgwas/quail


# format phenotype data
Rscript --vanilla vgwas/format_phenos_quail.R ${phenoFile} ${pheno} 


# Run Step 1
Rscript vgwas/QUAIL/Step1_QUAIL_rank_score.R \
--pheno ../data/vgwas/quail/ukb_${pheno}_quail.txt \
--covar ../data/vgwas/quail/ukb_covars_quail.txt \
--output ../data/vgwas/quail/ukb_${pheno}_rank_score.txt \
--num_levels 2000 \
--num_cores 5



# Run Step 2 as array
qsub -t 1-22 vgwas/vgwas.sh 0.005 ${pheno} 



#EOF
