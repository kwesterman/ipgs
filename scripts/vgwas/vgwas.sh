#!/bin/bash

#$ -l h_vmem=80G
#$ -l h_rt=12:00:00
#$ -o ../reports


#$ -j y
#$ -cwd



CHR=$SGE_TASK_ID

MAF=$1
pheno=$2


scratch=/broad/hptmp/gervis


source /broad/software/scripts/useuse
use R-4.1



# Run QUAIL Step 2 – vGWAS 
Rscript ../opt/QUAIL/Step2_QUAIL_vQTL_JEG.R \
--pheno ../data/processed/vgwas/${pheno}_pheno.txt \
--pheno_rs ../data/processed/vgwas/${pheno}_rankscore.txt \
--covar ../data/processed/vgwas/${pheno}_covars.txt \
--geno ${scratch}/plinkset/chr${CHR}_sel_maf${MAF} \
--output ../data/processed/vgwas/${pheno}_chr${CHR} \
--plink_path ../opt/plink2


## EOF
