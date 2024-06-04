#!/bin/bash


#$ -l h_vmem=10G
#$ -l h_rt=1:00:00

#$ -cwd
#$ -j y


tag=$1
pgs_type=$2


pgs_dir=../data/processed/pgs


source /broad/software/scripts/useuse
use R-4.1


# Load PGS weights and merge with clumped SNP set

R --vanilla <<EOF
library(tidyverse)
source("pgs/export_pgs_weights_funcs.R")
pgs_weights_df <- export_${pgs_type}("${tag}")
pgs_weights_df %>%
  mutate(CHR = paste0("chr", CHR)) %>%
  select(chrom = CHR, chromStart = POS, chromEnd = POS) %>%
  write_tsv("${pgs_dir}/${tag}_pgs_snps_hg19.txt", col_names = FALSE)
saveRDS(pgs_weights_df, "tmp_${tag}_pgs_weights.rds")
EOF


# Lift SNP positions from h19 to hg38

liftover=../opt/liftover/liftOver
${liftover} \
	${pgs_dir}/${tag}_pgs_snps_hg19.txt \
	../opt/liftover/hg19ToHg38.over.chain.gz \
	${pgs_dir}/${tag}_pgs_snps_hg38.txt \
	${pgs_dir}/${tag}_pgs_snps_hg38.unmapped

if [ -s "${pgs_dir}/${tag}_pgs_snps_hg38.unmapped" ]; then exit 1; fi


# Output final set of PGS weights

R --vanilla <<EOF
library(tidyverse)
pgs_weights_df <- readRDS("tmp_${tag}_pgs_weights.rds")
snp_df <- read_tsv("${pgs_dir}/${tag}_pgs_snps_hg38.txt",
		   col_names = c("CHR", "POS", "POS2"))
snp_df %>%
  select(CHR, POS) %>%
  mutate(CHR = gsub("chr", "", CHR)) %>%
  bind_cols(select(pgs_weights_df, -CHR, -POS)) %>%
  write_tsv("${pgs_dir}/${tag}_pgs_weights_hg38.txt")
EOF

rm tmp_${tag}_pgs_weights.rds
