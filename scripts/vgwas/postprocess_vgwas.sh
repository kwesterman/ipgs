#!/bin/bash


#$ -l h_vmem=15G
#$ -l h_rt=0:30:00
#$ -o ../reports/

#$ -cwd
#$ -j y


y=$1
fn_prefix=../data/processed/vgwas/${y}
MAF=0.005

ukb_bgen_dir=/broad/ukbb/imputed_v3 



# Merge chromosome-specific summary statistics
head -1 ${fn_prefix}_chr1_QUAIL_vQTL.txt > ${fn_prefix}_vQTL_merged
for chr in {1..22}; do
	echo "${fn_prefix}_chr${chr}..."
	tail -n +2 ${fn_prefix}_chr${chr}_QUAIL_vQTL.txt >> ${fn_prefix}_vQTL_merged
done


# Run postprocessing script
source /broad/software/scripts/useuse
use R-4.1

Rscript vgwas/postprocess_vgwas.R ${fn_prefix}_vQTL_merged ${y}

#EOF
