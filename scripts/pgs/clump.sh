#!/bin/sh


#$ -l os=RedHat7
#$ -l h_vmem=10G
#$ -l h_rt=1:00:00

#$ -j y
#$ -cwd


tag=$1

pgs_dir=../data/processed/pgs
input_file=${pgs_dir}/${tag}_pgsInput

ld_ref_prefix=../data/processed/ld_ref/ukb_20k_hg19

plink2=../opt/plink2

# Perform clumping

${plink2} \
	--pfile ${ld_ref_prefix} \
	--rm-dup force-first \
	--clump ${input_file} \
	--clump-p1 0.05 \
	--clump-r2 0.1 \
	--clump-kb 250 \
	--clump-id-field SNP \
	--clump-p-field P \
	--out ${pgs_dir}/${tag}

# Generate a file specifying P&T thresholds

pt_range_list_file=${pgs_dir}/${tag}_pt_range_list.txt
> ${pt_range_list_file}
for thresh in 0.05 0.01 0.00001 0.00000005; do 
        echo "${thresh} 0 ${thresh}" >> ${pt_range_list_file}; 
done

# List all SNPs that might be used in downstream PRS

awk 'NR!=1 {print $1}' ${input_file} > ${pgs_dir}/${tag}_pgs_snps.txt
