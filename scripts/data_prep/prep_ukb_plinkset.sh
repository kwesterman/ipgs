#!/bin/bash

#$ -l h_vmem=50G
#$ -l h_rt=2:00:00

#$ -cwd
#$ -j y


CHR=$SGE_TASK_ID
MAF=$1


scratch=/broad/hptmp/gervis
opt=../opt

ukb_bgen_dir=/broad/ukbb/imputed_v3 
ukb_sample_dir=/humgen/florezlab/UKBB_app27892 


source /broad/software/scripts/useuse

reuse -q Anaconda3
source activate $opt/bgen



# select snps

awk -v CHR=$CHR -v MAF=$MAF '{ if ($6>MAF && $8>0.3 ) {print $2} }' ${ukb_bgen_dir}/ukb_mfi_chr${CHR}_v3.txt > ${scratch}/snplist_chr${CHR}_maf${MAF}.txt

bgenix -g $ukb_bgen_dir/ukb_imp_chr${CHR}_v3.bgen -incl-rsids ${scratch}/snplist_chr${CHR}_maf${MAF}.txt > ${scratch}/chr${CHR}_sel_maf${MAF}.bgen



# bgen to bed/bim/fam

mkdir -p ${scratch}/plinkset

${opt}/plink2 --bgen ${scratch}/chr${CHR}_sel_maf${MAF}.bgen ref-first \
--sample ${ukb_sample_dir}/ukb27892_imp_chrAUT_v3_s487395.sample \
--make-bed \
--memory 50000 \
--rm-dup force-first \
--out ${scratch}/plinkset/chr${CHR}_sel_maf${MAF} \
&& rm ${scratch}/plinkset/chr${CHR}_sel_maf${MAF}.bgen \
&& rm ${scratch}/snplist_chr${CHR}_maf${MAF}.txt


#EOF


