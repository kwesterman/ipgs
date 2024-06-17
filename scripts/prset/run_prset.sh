#!/bin/sh


#$ -l os=RedHat7
#$ -l h_vmem=10G
#$ -l h_rt=6:00:00


#$ -pe smp 8
#$ -binding linear:8
#$ -R y

#$ -j y
#$ -cwd


tag=$1

pgs_dir=../data/processed/pgs
input_file=${pgs_dir}/${tag}_pgsInput
ld_ref_prefix=../data/processed/ld_ref/ukb_20k_hg19
plink2=../opt/plink2
prsice_dir=../opt/PRSice
prsice_datadir=../data/processed/prsice


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
        --bar-levels 5e-3 \
        --fastscore \
        --no-full \
        \
        --target /broad/ukbb/imputed_v3/ukb_imp_chr#_v3,/humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
        --type bgen \
	--exclude ${prsice_datadir}/ukb_duplicate_snpids.txt \
        --all-score \
        --no-regress \
        \
	--gtf ${prsice_datadir}/Homo_sapiens.GRCh37.75.gtf.gz \
	--msigdb ${prsice_datadir}/h.all.v2023.2.Hs.symbols.gmt \
	--set-perm 1 \
	--wind-3 1K \
	--wind-5 2K \
	\
	--print-snp \
	--out ${pgs_dir}/${tag}_prset \
	\
	--seed 123 \
        --thread 8 \
        --memory 10Gb


##Rscript ${prsice_dir}/PRSice.R --dir ${pgs_dir} \
##	--prsice ${prsice_dir}/PRSice_linux \
##	--base ${input_file} \
##	--A2 NEA \
##	--target ${ld_ref_prefix} \
##	--thread 1 \
##	--snp SNP \
##	--chr CHR \
##	--bp POS \
##	--A1 EA \
##	--stat beta \
##	--pvalue P \
##	--beta \
##	--binary-target F \
##	--no-full \
##	--lower 5e-8 \
##	--upper 0.05 \
##	--interval 5e-2 \
##	--no-regress \
##	--gtf Homo_sapiens.GRCh37.75.gtf.gz \
##	--msigdb ${prset_datadir}/h.all.v2023.2.Hs.symbols.gmt \
##	--set-perm 1 \
##	--wind-3 1K \
##	--wind-5 2K
