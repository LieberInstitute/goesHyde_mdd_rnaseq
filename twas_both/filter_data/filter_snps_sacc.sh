#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -pe local 4
#$ -N "filter_snps_MDD_genes_sacc"
#$ -j y
#$ -o logs/filter_snps_sacc_gene_$JOB_ID.txt

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load dependencies
module load plink/1.90b6.6
module load conda_R/4.0.x

## List current modules
module list

mkdir logs

Rscript filter_snps.R -r "sACC"
