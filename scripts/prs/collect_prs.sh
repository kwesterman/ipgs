#!/bin/bash


#$ -l h_vmem=50G
#$ -l h_rt=2:00:00

#$ -cwd
#$ -j y

tag=$1

source /broad/software/scripts/useuse
use R-4.1

Rscript prs/collect_prs.R ${tag}
