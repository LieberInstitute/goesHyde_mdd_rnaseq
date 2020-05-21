#!/bin/bash
#$ -cwd
#$ -l mem_free=80G,h_vmem=80G,h_fsize=100G
#$ -N qSV_model_DE_analysis
#$ -o logs/qSV_model_DE_analysis.txt
#$ -e logs/qSV_model_DE_analysis.txt
#$ -hold_jid get_expression_cutoffs,merge_MDDseqBiPseq_degradation
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
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
Rscript qSV_model_DE_analysis.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/