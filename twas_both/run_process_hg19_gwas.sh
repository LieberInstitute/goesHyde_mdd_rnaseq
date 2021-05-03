#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -pe local 1
#$ -N process-hg19-gwas
#$ -j y
#$ -o logs/hg19_gwas.txt

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load dependencies
module load conda_R/4.0.x

## List current modules
module list

Rscript process-hg19-gwas.R

echo "**** Job ends ****"
date
