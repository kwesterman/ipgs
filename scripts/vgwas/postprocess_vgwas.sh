#!/bin/bash


#$ -l h_vmem=15G
#$ -l h_rt=0:30:00

#$ -cwd
#$ -j y


y=$1
fn_prefix=../data/processed/vgwas/${y}
ukb_bgen_dir=/broad/ukbb/imputed_v3 



# Merge chromosome-specific summary statistics & create list of high quality variants
head -1 ${fn_prefix}_chr1 > ${fn_prefix}_merged
for chr in {1..22}; do
	echo "${fn_prefix}_chr${chr}..."
	tail -n +2 ${fn_prefix}_chr${chr} >> ${fn_prefix}_merged
done


# Create list of high quality variants (info >0.5)
echo "" > ../data/processed/ukb_geno_1to22_maf${MAF}_info0.5.txt
for chr in {1..22}; do
	echo "Filtering chromosome ${chr}..."
	awk -v chr=${chr} -v maf=${MAF} '$6 > maf && $8 > 0.5 {print chr"\t"$0}' ${ukb_bgen_dir}/ukb_mfi_chr${chr}_v3.txt >> ../data/processed/ukb_geno_1to22_maf${MAF}_info0.5.txt
done


# Run postprocessing script
source /broad/software/scripts/useuse
use R-4.1

Rscript vgwas/postprocess_vgwas.R ${fn_prefix}_merged ${y}

#EOF
11826013