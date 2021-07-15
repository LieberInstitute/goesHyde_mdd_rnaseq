#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -N run_PEER
#$ -j y
#$ -o /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/logs/run_PEER.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${TASK_ID}"

cd "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/"

conda activate r-peer

num_peer=60

prefix_Amygdala="goesHyde_mdd_rnaseq_Amygdala"

echo "Begin PEER script"
time {
  Rscript code/01_get_inv_quantile_norm/run_PEER.R processed-data/01_get_inv_quantile_norm/${prefix_Amygdala}.expression.bed.gz ${prefix_Amygdala} ${num_peer}
}
echo "PEER script ended"

mv goesHyde_mdd_rnaseq* processed-data/01_get_inv_quantile_norm/

echo "**** Job ends ****"
date
