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

# create junc_count_files.txt with input_check.R

cat ../data/junc_count_files.txt | while read line; 
	do python reformat_junc.py $line;
done

ls ../data/junc/*.junc > ../data/all_junc_filenames.txt

echo "**** Job ends ****"
date
