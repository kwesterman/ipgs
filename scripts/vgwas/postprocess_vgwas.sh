#!/bin/bash

#$ -l h_vmem=4G
#$ -l h_rt=0:10:00

#$ -j y
#$ -cwd


pheno=$1

source /broad/software/scripts/useuse
use R-4.1


R --vanilla <<EOF
library(tidyverse) ; library(data.table)
do.call(rbind.data.frame, lapply(as.list(1:22), function(i) {fread(paste0("../data/vgwas/quail/ukb_chr", i, "_${pheno}_rank_score_QUAIL_vQTL.txt")) } )) %>% fwrite(paste0("../data/vgwas/ukb_${pheno}_QUAIL_vQTL.txt"))
EOF


#EOF
