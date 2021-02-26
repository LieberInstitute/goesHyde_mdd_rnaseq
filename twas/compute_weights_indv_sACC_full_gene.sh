#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=4G,h_vmem=4G,h_fsize=100G
#$ -N compute_weights_indv_sACC_full_gene
#$ -j y
#$ -o logs/sACC_gene/compute_weights_indv_sACC_full_gene.$TASK_ID.txt
#$ -t 1-24958
#$ -tc 40
#$ -m a

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load dependencies
module load plink/1.90b6.6
module load fusion_twas/github
module load conda_R/4.0.x

## List current modules
module list

# relative path for FILELIST
FILELIST=$(echo "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/twas/sACC_gene/input_ids.txt")

## File id and feature name
FEATURENUM=$(awk 'BEGIN {FS="\t"} {print $1}' ${FILELIST} | awk "NR==${SGE_TASK_ID}")
FEATUREID=$(awk 'BEGIN {FS="\t"} {print $2}' ${FILELIST} | awk "NR==${SGE_TASK_ID}")

## Define files
FILTBIM="/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/twas/sACC_gene/bim_files/sACC_gene_${FEATURENUM}/filtered_snps_sACC_gene_${FEATURENUM}"
TMPFILES="tmp_files/gene_${FEATURENUM}"
OUTFILES="out_files/gene_${FEATURENUM}"

cd sACC_gene/

Rscript /jhpce/shared/jhpce/libd/fusion_twas/github/fusion_twas/FUSION.compute_weights.R \
    --bfile ${FILTBIM} \
    --tmp ${TMPFILES} \
    --out ${OUTFILES} \
    --PATH_gemma /dcl01/lieber/ajaffe/lab/twas/software/gemma-0.98.1-linux-static \
    --models top1,blup,lasso,enet --hsq_p 1.0001 --verbose 1 --save_hsq

echo "**** Job ends ****"
date
