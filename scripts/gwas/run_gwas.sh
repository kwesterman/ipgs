#!/bin/sh


#$ -l os=RedHat7
#$ -l h_vmem=5G
#$ -l h_rt=12:00:00

#$ -pe smp 8
#$ -binding linear:8
#$ -R y

#$ -j y
#$ -cwd


pheno=$1

chr=$SGE_TASK_ID


gPC_arr=(gPC1 gPC2 gPC3 gPC4 gPC5 gPC6 gPC7 gPC8 gPC9 gPC10)
covars="sex age age_squared ageBySex ${gPC_arr[@]}"


source /broad/software/scripts/useuse
use R-4.1
R --no-save <<EOF
library(tidyverse)
read_csv("../data/processed/ukb_training_set.csv") %>%
  select(where(~ !is.character(.x))) %>%
  write_csv("../data/processed/${pheno}_phenos_chr${chr}.tmp")
EOF

singularity_dir=~/kw/singularity
singularity exec \
	-B ../data/processed:/data \
	-B /broad/ukbb/imputed_v3:/bgendir \
	-B /humgen/florezlab/UKBB_app27892:/sampledir \
	${singularity_dir}/gem-v1.5.2-workflow.simg \
	/bin/bash <<EOF

/GEM/GEM \
	--bgen /bgendir/ukb_imp_chr${chr}_v3.bgen \
	--sample /sampledir/ukb27892_imp_chrAUT_v3_s487395.sample \
	--pheno-file /data/${pheno}_phenos_chr${chr}.tmp \
	--sampleid-name id \
	--pheno-name ${pheno} \
	--covar-names ${covars} \
	--delim , \
	--missing-value NA \
	--cat-threshold 3 \
	--maf 0.01 \
	--robust 1 \
	--threads 8 \
	--output-style meta \
	--out /data/gwas/${pheno}_chr${chr}

EOF

rm ../data/processed/${pheno}_phenos_chr${chr}.tmp
