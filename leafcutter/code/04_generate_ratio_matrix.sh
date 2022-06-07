#!/bin/bash
#$ -cwd
#$ -l mem_free=25G,h_vmem=25G,h_fsize=100G
#$ -N generate_ratio_matrix
#$ -o logs/generate_ratio_matrix.txt
#$ -e logs/generate_ratio_matrix.txt

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

module load python/2.7

# Create qqnormed file
echo "build qqnorm file..."
python prepare_phenotype_table.py  ../data/clusters/leafcutter_perind.counts.gz -p 10

echo "**** Job ends ****"
date