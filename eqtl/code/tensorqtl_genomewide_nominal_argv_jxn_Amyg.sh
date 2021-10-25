#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -N tensorqtl_genomewide_nominal_argv_jxn_Amyg
#$ -o logs/tensorqtl_genomewide_nominal_argv_jxn_Amyg.txt
#$ -e logs/tensorqtl_genomewide_nominal_argv_jxn_Amyg.txt
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
conda activate hello
## get pair

python tensorqtl_genomewide_nominal_argv.py jxn_Amygdala

echo "**** Job ends ****"
date
