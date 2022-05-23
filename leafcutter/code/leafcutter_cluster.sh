#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N leafcutter_cluster
#$ -o logs/leafcutter_cluster.txt
#$ -e logs/leafcutter_cluster.txt

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

python leafcutter_cluster.py -j ../data/all_jxn_filenames.txt -m 50 -o BP_RNAseq -l 500000

echo "**** Job ends ****"
date