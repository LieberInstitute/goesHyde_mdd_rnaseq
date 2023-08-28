#!/bin/bash

## Usage:
# sh 03_tensorqtl_genomewide_nominal.sh

## Create the logs directory
mkdir -p logs

for feature_region in gene_Amygdala gene_sACC exon_Amygdala exon_sACC jxn_Amygdala jxn_sACC tx_Amygdala tx_sACC; do

    ## Internal script name
    SHORT="03_tensorqtl_genomewide_nominal_${feature_region}"
    NAME="tensorqtl_genomewide_nominal_${feature_region}"

    # Construct shell file
    echo "Creating script 03_tensorqtl_genomewide_nominal_${feature_region}"
    cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l caracol,mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -N ${NAME}
#$ -o logs/${SHORT}.txt
#$ -e logs/${SHORT}.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
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

## Run Python
python 03_tensorqtl_genomewide_nominal.py ${feature_region}
echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


EOF

    call="qsub .${SHORT}.sh"
    echo $call
    $call
done
