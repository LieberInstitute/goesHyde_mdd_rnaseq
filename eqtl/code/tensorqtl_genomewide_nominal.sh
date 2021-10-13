#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -N tensorqtl_genomewide_nominal
#$ -o logs/tensorqtl_genomewide_nominal.txt
#$ -e logs/tensorqtl_genomewide_nominal.txt
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

## List current modules for reproducibility
module list

## Edit with your job command
conda activate hello
echo "Hello!"
python tensorqtl_genomewide_nominal.py

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
