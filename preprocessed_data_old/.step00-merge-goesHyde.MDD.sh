#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N step00-merge
#$ -pe local 20
#$ -o logs/merge-goesHyde.txt
#$ -e logs/merge-goesHyde.txt
#$ -m a

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

Rscript /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/sh/step00-merge.R -s /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/.samples.manifest_unmerged -o /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/FASTQ_merged -c 24

echo "**** Job ends ****"
date
