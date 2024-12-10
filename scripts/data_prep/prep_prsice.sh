#!/bin/bash


#$ -l h_vmem=10G

#$ -cwd
#$ -j y


prsice_dir=../data/processed/prsice
ukb_imp_dir=/broad/ukbb/imputed_v3
ldref_dir=../data/processed/ld_ref

source /broad/software/scripts/useuse
use R-4.1


# Install PRSice

mkdir -p ../opt/PRSice \
	&& wget https://github.com/choishingwan/PRSice/releases/download/2.3.5/PRSice_linux.zip \
	&& mv PRSice_linux.zip ../opt/PRSice/ \
	&& unzip ../opt/PRSice/PRSice_linux.zip

# Generate list of duplicated UKB SNP IDs

"" > ${prsice_dir}/ukb_duplicate_snpids.txt
for chr in {1..22}; do
	echo "Finding duplicate rsIDs in chromosome ${chr}..."
	cat ${ukb_imp_dir}/ukb_mfi_chr${chr}_v3.txt | cut -f2 | sort | uniq -d >> ${prsice_dir}/ukb_duplicate_snpids.txt
done
echo "Finding duplicate rsIDs in LD reference file..."
cat ${ldref_dir}/ukb_20k_hg19.bim | cut -f2 | sort | uniq -d >> ${prsice_dir}/ukb_duplicate_snpids.txt
