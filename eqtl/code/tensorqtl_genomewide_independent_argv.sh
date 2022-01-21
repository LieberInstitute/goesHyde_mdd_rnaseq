#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -pe local 1
#$ -N tensorqtl_genomewide_nominal_argv
#$ -o logs/tensorqtl_genomewide_independent_argv.$TASK_ID.txt
#$ -e logs/tensorqtl_genomewide_independent_argv.$TASK_ID.txt
#$ -t 1-46
#$ -tc 4

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

## get pair
PAIR=$(awk "NR==${SGE_TASK_ID}" region_chr.txt)
echo "Processing pair: ${PAIR}"

python tensorqtl_genomewide_independent_argv.py ${PAIR}

echo "**** Job ends ****"
date
