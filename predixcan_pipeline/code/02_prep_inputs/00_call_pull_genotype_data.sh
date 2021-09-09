#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=200G,h_vmem=200G,h_fsize=200G
#$ -N pull_genotype_data
#$ -j y
#$ -o /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/logs/pull_genotype_data.log

# This is a bespoke qsub shell file written simply to call the 02_pull_genotype_data Rscript
# NOT AN OFFICIAL PART OF THE PIPELINE at least not yet

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

Rscript /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/code/02_prep_inputs/01_pull_genotype_data.R

echo "**** Job ends ****"
date
