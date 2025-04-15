#!/bin/bash


#$ -l h_vmem=10G
#$ -l h_rt=6:00:00

#$ -cwd
#$ -j y


source /broad/software/scripts/useuse
use GCC-5.2
use R-4.0


sumstats_file=$1
prefix=$2


magma_dir=../opt/magma_v1.10
data_dir=../data/raw/magma
ldref_dir=../data/processed/ld_ref
annot_dir=../data/processed/magma
output_dir=../data/processed/magma

geneloc_file=${data_dir}/NCBI37.3.gene.loc


# Gene analysis (based on p-values)
${magma_dir}/magma \
	--bfile ${ldref_dir}/ukb_20k_hg19 \
	--pval ${sumstats_file} use=SNP,P_marg N=350000 \
	--gene-model snp-wise=mean \
	--gene-annot ${annot_dir}/ukb_20k_hg19_2.1.genes.annot \
	--out ${output_dir}/${prefix}

pathway_file=${annot_dir}/c2.all.v2024.1.Hs.entrez.gmt
${magma_dir}/magma \
    --gene-results ${output_dir}/${prefix}.genes.raw \
    --set-annot ${pathway_file} \
    --out ${output_dir}/${prefix}
