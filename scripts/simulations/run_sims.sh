#!/bin/bash


#$ -l h_vmem=3G
#$ -l h_rt=12:00:00

#$ -pe smp 16
#$ -binding linear:16
#$ -R y

#$ -cwd
#$ -j y

source /broad/software/scripts/useuse
use R-4.1

Rscript simulations/run_sims.R
