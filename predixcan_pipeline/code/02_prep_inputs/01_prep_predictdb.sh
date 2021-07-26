#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -N PredictDB_prep
#$ -j y
#$ -o /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/logs/PredictDB_prep.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

conda activate eqtl_prepare_expression

python /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/PredictDB-Tutorial/code/parse_gtf.py \
  /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf \
  /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/02_prep_inputs/gencode.v25.annotationGRCh38_PARSED.gtf

echo "**** Job ends ****"
date
