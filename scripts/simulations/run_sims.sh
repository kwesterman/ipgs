#!/bin/bash


#$ -l h_vmem=3G
#$ -l h_rt=12:00:00

#$ -pe smp 32
#$ -binding linear:32
#$ -R y

#$ -cwd
#$ -j y

source /broad/software/scripts/useuse
use R-4.1

Rscript simulations/run_sims.R 32
