#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -N extract_qqnorm_data
#$ -o logs/extract_qqnorm_data.txt
#$ -e logs/extract_qqnorm_data.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R

## List current modules for reproducibility
module list

# Create file for Amygdala and sACC regions containing qqnorm'd data with Brain Numbers (BrNum) as column names
Rscript extract_qqnorm_data.R

# Loading data into R changes the '#' in the header line to 'X.', the following sed commands revert the change
sed -i 's/X./#/' ../data/qqnorm/amygdala_qqnorm.txt
sed -i 's/X./#/' ./data/qqnorm/sacc_qqnorm.txt

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
