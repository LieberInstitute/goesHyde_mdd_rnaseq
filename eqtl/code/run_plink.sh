#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N run_plink
#$ -o logs/run_plink.txt
#$ -e logs/run_plink.txt

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load plink

## List current modules for reproducibility
module list

## Edit with your job command
plink --vcf LIBD_maf01_gwas_BPD.vcf.gz --make-bed --out LIBD_maf01_gwas_BPD

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
