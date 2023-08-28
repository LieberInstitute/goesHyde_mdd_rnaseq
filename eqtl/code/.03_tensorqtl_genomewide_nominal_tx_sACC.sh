#!/bin/bash
#$ -cwd
#$ -l caracol,mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -N tensorqtl_genomewide_nominal_tx_sACC
#$ -o logs/03_tensorqtl_genomewide_nominal_tx_sACC.txt
#$ -e logs/03_tensorqtl_genomewide_nominal_tx_sACC.txt
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
module load tensorqtl
module list

USAGE_CUTOFF=10
NUM_GPUS=1

avail_gpus=

#  Simply exit with an error if there are no GPUs left
if [[ -z  ]]; then
    echo "No GPUs are available."
    exit 1
fi

export CUDA_VISIBLE_DEVICES=

## Run Python
python 03_tensorqtl_genomewide_nominal.py tx_sACC
echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


