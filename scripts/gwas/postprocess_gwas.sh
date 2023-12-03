#!/bin/bash


#$ -l h_vmem=15G
#$ -l h_rt=00:30:00

#$ -cwd
#$ -j y


y=$1
fn_prefix=../data/processed/gwas/${y}


# Merge chromosome-specific summary statistics
head -1 ${fn_prefix}_chr1 > ${fn_prefix}_merged
for chr in {1..22}; do
	echo "${fn_prefix}_chr${chr}..."
	tail -n +2 ${fn_prefix}_chr${chr} >> ${fn_prefix}_merged
done


# Run postprocessing script
source /broad/software/scripts/useuse
use R-4.1

Rscript gwas/postprocess_gwas.R ${fn_prefix}_merged ${y}
