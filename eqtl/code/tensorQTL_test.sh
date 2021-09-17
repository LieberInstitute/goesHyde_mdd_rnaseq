#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N tensorQTL_test
#$ -o logs/tensorQTL_test.txt
#$ -e logs/tensorQTL_test.txt
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
## module load python

## List current modules for reproducibility
module list

## Edit with your job command
python3 -m tensorqtl ../data/risk_snps/LIBD_maf01_gwas_BPD ../data/phenotype_bed/gene_Amygdala.bed.gz ../data/tensorQTL_out/gene_Amygdala --covariates ../data/covariates_txt/covariates_gene_Amygdala.txt --mode cis

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
