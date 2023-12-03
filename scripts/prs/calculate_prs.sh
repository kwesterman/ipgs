#!/bin/sh


#$ -l os=RedHat7
#$ -l h_vmem=10G
#$ -l h_rt=1:00:00

###$ -pe smp 8
###$ -binding linear:8
###$ -R y

#$ -j y
#$ -cwd


tag=$1

chr=$SGE_TASK_ID

prs_dir=../data/processed/prs
input_file=${prs_dir}/${tag}_prsInput

ld_ref_prefix=../data/processed/ld_ref/ukb_20k_hg19

source /broad/software/scripts/useuse
use Anaconda3 

plink2=../opt/plink2

# Extract a plinkset with all potential PRS variants

#awk 'NR!=1 {print $1}' ${input_file} > ${prs_dir}/${tag}_prs_snps.txt

source activate bgen
bgenix \
	-g /broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen \
	-incl-rsids ${prs_dir}/${tag}_prs_snps.txt \
	> ${prs_dir}/${tag}_prs_subset_chr${chr}.bgen

${plink2} \
	--bgen ${prs_dir}/${tag}_prs_subset_chr${chr}.bgen ref-first \
	--sample /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
	--extract ${prs_dir}/${tag}_prs_snps.txt \
	--rm-dup force-first \
	--make-pgen \
	--out ${prs_dir}/${tag}_prs_subset_chr${chr} \
	--memory 5000

##${plink2} \
##        --bgen /broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen ref-first \
##        --sample /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
##        --extract ${prs_dir}/${tag}_prs_snps.txt \
##	--rm-dup force-first \
##        --make-pgen \
##	--out ${prs_dir}/${tag}_prs_subset_chr${chr} \
##        --memory 40000 \
##	--threads 8

# Calculate a series of P&T-based scores

${plink2} \
	--pfile ${prs_dir}/${tag}_prs_subset_chr${chr} \
        --score ${input_file} 1 4 6 header ignore-dup-ids cols=scoresums\
        --q-score-range ${prs_dir}/${tag}_pt_range_list.txt ${input_file} 1 8 header \
	--out ${prs_dir}/${tag}_chr${chr} \
        --memory 5000
