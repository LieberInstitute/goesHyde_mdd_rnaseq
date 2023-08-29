#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -N tensorqtl_genomewide_cis
#$ -o logs/05_tensorqtl_genomewide_cis.txt
#$ -e logs/05_tensorqtl_genomewide_cis.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load tensorqtl

## List current modules for reproducibility
module list

## Edit with your job command
python 05_tensorqtl_genomewide_cis.py

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/