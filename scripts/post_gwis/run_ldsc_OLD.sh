#!/bin/sh


#$ -l h_vmem=10G
#$ -l h_rt=1:00:00
#$ -cwd
#$ -j y


source /broad/software/scripts/useuse
use Anaconda3

source activate ldsc


ldsc_dir=../data/processed/ldsc
ldref_dir=../data/processed/ld_ref


tag1=$1
tag2=$2

working_dir=${tag1}_${tag2}_ldsc_dir


# Munge summary statistics to prep for LDSC
mkdir -p ${working_dir}
for tag in ${tag1} ${tag2}; do 
cat ${ldsc_dir}/${tag}_ldscInput > ${working_dir}/${tag}_ss
../opt/ldsc/munge_sumstats.py \
        --sumstats ${working_dir}/${tag}_ss \
        --merge-alleles ../data/raw/ldsc/w_hm3.snplist \
        --snp SNP \
        --N-col N \
        --a1 EA \
        --a2 NEA \
        --p P \
        --frq AF \
        --signed-sumstats Beta,0 \
        --chunksize 500000 \
        --out ${working_dir}/${tag}
rm ${working_dir}/${tag}_ss
done


# Run LDSC to calculate genetic correlations
../opt/ldsc/ldsc.py \
	--rg ${working_dir}/${tag1}.sumstats.gz,${working_dir}/${tag2}.sumstats.gz \
	--ref-ld ${ldref_dir}/ukb_20k_hg19_withCM \
	--w-ld ${ldref_dir}/ukb_20k_hg19_withCM \
	--out ${ldsc_dir}/${tag1}_${tag2}

rm -r ${working_dir}

conda deactivate
