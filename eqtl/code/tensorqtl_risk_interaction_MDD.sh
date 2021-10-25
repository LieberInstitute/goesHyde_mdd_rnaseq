#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -N tensorqtl_risk_interaction
#$ -o logs/tensorqtl_risk_interaction_MDD.txt
#$ -e logs/tensorqtl_risk_interaction_MDD.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

conda activate hello
## List current modules for reproducibility
module list

## Edit with your job command
python tensorqtl_risk_interaction_MDD.py

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
