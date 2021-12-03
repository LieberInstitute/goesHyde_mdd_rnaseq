#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -N variancePartition_exploration
#$ -o logs/variancePartition_exploration.$TASK_ID.txt
#$ -e logs/variancePartition_exploration.$TASK_ID.txt
#$ -m e
#$ -t 1-4
#$ -tc 4

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

## get pair
PAIR=$(awk "NR==${SGE_TASK_ID}" region_dx.txt)
echo "Processing pair: ${PAIR}"

## Edit with your job command
Rscript variancePartition_exploration.R ${PAIR}

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
