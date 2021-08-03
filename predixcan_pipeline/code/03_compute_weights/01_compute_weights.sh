#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -N PredictDB_compute_weights
#$ -j y
#$ -o /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/logs/PredictDB_compute_weights.$TASK_ID.log
#$ -t 1-22

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

module load conda_R/4.1

echo "Computing weights for chromosome $TASK_ID"
time Rscript /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/PredictDB-Tutorial/code/gtex_tiss_chrom_training.R ${TASK_ID}

echo "**** Job ends ****"
date
