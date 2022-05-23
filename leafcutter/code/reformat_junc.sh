#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N reformat_junc
#$ -o logs/reformat_junc.txt
#$ -e logs/reformat_junc.txt

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

for i in $(ls /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/Counts/junction/*.count); 
	do python reformat_junc.py $i;
done

ls ../data/jucn/*.junc >> ../data/all_jxn_filenames.txt

echo "**** Job ends ****"
date