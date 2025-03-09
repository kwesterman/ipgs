#!/bin/bash

#$ -l h_vmem=30G
#$ -l h_rt=6:00:00
#$ -o ../reports

#$ -j y
#$ -cwd


pheno=$1
MAF=0.005
phenoFile=/humgen/florezlab/UKBB_app27892/UKBB_app27892_kew/ipgs/data/processed/ukb_training_set.csv

#phenoFile=../data/processed/ukb_training_set_TEST.csv

source /broad/software/scripts/useuse
use R-4.1


mkdir -p ../reports



# prep pheno & covar file for quail
Rscript --vanilla vgwas/preprocess_vgwas.R ${phenoFile} ${pheno}


# use quail to calculate integrated rank score
Rscript ../opt/QUAIL/Step1_QUAIL_rank_score.R \
--pheno ../data/processed/vgwas/${pheno}_pheno.txt \
--covar ../data/processed/vgwas/${pheno}_covars.txt \
--output ../data/processed/vgwas/${pheno}_rankscore.txt \
--num_levels 2000 \
--num_cores 5



# Run Step 2 as array
use UGER
qsub -t 1-22 vgwas/vgwas.sh ${pheno}


#EOF

