#!/bin/bash

#$ -l h_vmem=8G
#$ -l h_rt=1:00:00

#$ -j y
#$ -cwd

MAF=$1
scratch=/broad/hptmp/gervis


# delete bgen files
rm ${scratch}/chr*_sel_maf${MAF}.bgen


##EOF
