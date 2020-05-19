#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=3G,h_vmem=3G,h_fsize=100G
#$ -pe local 8
#$ -N run_wgcna_combined
#$ -o logs/run_wgcna_combined.txt
#$ -e logs/run_wgcna_combined.txt
#$ -hold_jid get_expression_cuttoff,merge_MDDseq_BiPseq_degradation
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
module load conda_R/3.6.x

## List current modules for reproducibility
module list

## Edit with your job command EDIT EDIT
Rscript run_wgcna_combined.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
