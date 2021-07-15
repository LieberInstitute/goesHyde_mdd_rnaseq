#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N combine_covs
#$ -j y
#$ -o /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/logs/combine_covs.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${TASK_ID}"

conda activate eqtl_prepare_expression

prefix_Amygdala="goesHyde_mdd_rnaseq_Amygdala"
genotype_pcs="/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/goesHyde_bipolarMdd_Genotypes_mds.csv"

Rscript code/01_get_inv_quantile_norm/03_convert_rda.R

./code/01_get_inv_quantile_norm/combine_covariates.py processed-data/01_get_inv_quantile_norm/${prefix_Amygdala}.PEER_covariates.txt ${prefix_Amygdala} \
    --genotype_pcs ${genotype_pcs}

echo "**** Job ends ****"
date
