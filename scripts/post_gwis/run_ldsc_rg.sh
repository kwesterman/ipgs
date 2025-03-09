#!/bin/sh


#$ -l h_vmem=10G
#$ -l h_rt=00:30:00
#$ -cwd
#$ -j y


bm=$1
ancestry=EUR


source /broad/software/scripts/useuse
use Anaconda3
use R-4.1


ss_dir=../data/processed/gwis
mkdir -p ${ss_dir}/ldsc	

source activate ldsc


# Download necessary reference LD files
if [ ! -e ../data/raw/ldsc/eur_w_ld_chr ]; then
	wget -P ../data/raw/ldsc/ https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
	bzip2 -d ../data/raw/ldsc/eur_w_ld_chr.tar.bz2
	tar xvf ../data/raw/ldsc/eur_w_ld_chr.tar -C ../data/raw/ldsc/ 
fi
if [ ! -f ../data/raw/ldsc/w_hm3.snplist ]; then
	wget -P ../data/raw/ldsc/ https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
	bzip2 -d ../data/raw/ldsc/w_hm3.snplist.bz2
fi

# Associate vQTL summary stats files with rsIDs
R --vanilla <<EOF
  library(tidyverse)
  sumstats <- data.table::fread("${ss_dir}/bmi_${bm}_merged", data.table = FALSE, stringsAsFactors = FALSE) %>%
    select(RSID, Effect_Allele, Non_Effect_Allele, AF, N_Samples, beta_main = Beta_G, SE_main = robust_SE_Beta_G, beta_int = all_of("Beta_G-bmi"), SE_int = all_of("robust_SE_Beta_G-bmi")) %>%
    mutate(P_main = pchisq(beta_main^2 / SE_main^2, df = 1, lower.tail = FALSE),
           P_int = pchisq(beta_int^2 / SE_int^2, df = 1, lower.tail = FALSE)) %>%
    select(RSID, Effect_Allele, Non_Effect_Allele, AF, N_Samples, beta_main, P_main, beta_int, P_int)
  sumstats %>%
    write_tsv("${ss_dir}/ldsc/bmi_${bm}_sumstats")
EOF

# Munge sumstats for main effect
../opt/ldsc/munge_sumstats.py \
	--sumstats ${ss_dir}/ldsc/bmi_${bm}_sumstats \
	--merge-alleles ../data/raw/ldsc/w_hm3.snplist \
	--snp RSID \
	--N-col N_Samples \
	--a1 Effect_Allele \
	--a2 Non_Effect_Allele \
	--p P_main \
	--frq AF \
	--signed-sumstats beta_main,0 \
	--chunksize 500000 \
	--out ${ss_dir}/ldsc/${bm}_main_ldsc

# Munge sumstats for interaction effect
../opt/ldsc/munge_sumstats.py \
	--sumstats ${ss_dir}/ldsc/bmi_${bm}_sumstats \
	--merge-alleles ../data/raw/ldsc/w_hm3.snplist \
	--snp RSID \
	--N-col N_Samples \
	--a1 Effect_Allele \
	--a2 Non_Effect_Allele \
	--p P_int \
	--frq AF \
	--signed-sumstats beta_int,0 \
	--chunksize 500000 \
	--out ${ss_dir}/ldsc/${bm}_int_ldsc

../opt/ldsc/ldsc.py \
	--rg ${ss_dir}/ldsc/${bm}_main_ldsc.sumstats.gz,${ss_dir}/ldsc/${bm}_int_ldsc.sumstats.gz \
	--ref-ld-chr ../data/raw/ldsc/eur_w_ld_chr/ \
	--w-ld-chr ../data/raw/ldsc/eur_w_ld_chr/ \
	--out ${ss_dir}/ldsc/${bm}

conda deactivate
