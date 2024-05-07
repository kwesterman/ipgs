#!/bin/bash

#$ -l h_vmem=80G
#$ -l h_rt=6:00:00
#$ -o ../reports

#$ -j y
#$ -cwd



MAF=$1
phenoFile=$2
pheno=$3


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
qsub -t 1-22 vgwas/vgwas.sh ${MAF} ${pheno} 



#EOF

