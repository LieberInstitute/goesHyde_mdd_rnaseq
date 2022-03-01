#!/bin/bash
#$ -cwd
#$ -l caracol,mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -N tensorqtl_genomewide_independent
#$ -o logs/tensorqtl_genomewide_independent.txt
#$ -e logs/tensorqtl_genomewide_independent.txt

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## List current modules for reproducibility
module load tensorqtl
module list

USAGE_CUTOFF=10
NUM_GPUS=1

avail_gpus=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader | cut -d " " -f 1 | awk -v usage="$USAGE_CUTOFF" '$1 < usage {print NR - 1}')

#  Simply exit with an error if there are no GPUs left
if [[ -z $avail_gpus ]]; then
    echo "No GPUs are available."
    exit 1
fi

export CUDA_VISIBLE_DEVICES=$(echo "$avail_gpus" | head -n $NUM_GPUS | paste -sd ",")

python tensorqtl_genomewide_independent.py ${PAIR}

echo "**** Job ends ****"
date
