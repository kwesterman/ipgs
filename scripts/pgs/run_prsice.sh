#!/bin/sh


#$ -l os=RedHat7
#$ -l h_vmem=10G
#$ -l h_rt=8:00:00


#$ -pe smp 8
#$ -binding linear:8
#$ -R y

#$ -j y
#$ -cwd


tag=$1

pgs_dir=../data/processed/pgs
input_file=${pgs_dir}/${tag}_pgsInput
ld_ref_prefix=../data/processed/ld_ref/ukb_20k_hg19
prsice_datadir=../data/processed/prsice
plink2=../opt/plink2
prsice_dir=../opt/PRSice


source /broad/software/scripts/useuse
use R-4.1


Rscript ${prsice_dir}/PRSice.R --dir ${pgs_dir} \
        --prsice ${prsice_dir}/PRSice_linux \
        \
        --base ${input_file} \
        --snp SNP \
        --chr CHR \
        --bp POS \
        --A1 EA \
        --A2 NEA \
        --stat beta \
        --pvalue P \
        --beta \
        \
        --ld ${ld_ref_prefix} \
	--clump-kb 250 \
	--clump-r2 0.1 \
        --bar-levels 5e-2,1e-2,5e-3,1e-3,5e-4,1e-4,5e-5,1e-5,5e-6,1e-6,5e-7,1e-7,5e-8 \
        --fastscore \
        --no-full \
        \
        --target /broad/ukbb/imputed_v3/ukb_imp_chr#_v3,/humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
        --type bgen \
	--exclude ${prsice_datadir}/ukb_duplicate_snpids.txt \
        --all-score \
        --no-regress \
        \
	--print-snp \
	--out ${pgs_dir}/${tag} \
	\
	--seed 123 \
        --thread 8 \
        --memory 10Gb
