#!/bin/bash

#$ -l h_vmem=80G
#$ -l h_rt=20:00:00
#$ -o ../reports


#$ -j y
#$ -cwd



CHR=$SGE_TASK_ID

pheno=$1


MAF=0.005
vgwas_dir=../data/processed/vgwas
scratch=/broad/hptmp/gervis

fn_prefix=${vgwas_dir}/${pheno}



source /broad/software/scripts/useuse
use R-4.1



# Run QUAIL Step 2 – vGWAS 
Rscript ../opt/QUAIL/Step2_QUAIL_vQTL_JEG.R \
--pheno ${fn_prefix}_pheno.txt \
--pheno_rs ${fn_prefix}_rankscore.txt \
--covar ${fn_prefix}_covars.txt \
--geno ${scratch}/plinkset/chr${CHR}_sel_maf${MAF} \
--output ${fn_prefix}_chr${CHR} \
--plink_path ../opt/plink2


## EOF
