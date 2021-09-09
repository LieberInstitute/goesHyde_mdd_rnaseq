#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=100G
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

module load conda_R/4.1

# Rscript /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/code/02_prep_inputs/02_pull_genotype_data.R

mkdir /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/02_prep_inputs/split_geno

Rscript /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/code/02_prep_inputs/03_process_snp_anno.R

# removing the redundant chromosome column
# for genos in ls /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/02_prep_inputs/split_geno/*; do
#   echo $genos;
#   awk '{NF=""; print $0}' $genos > $genos
# done

conda activate eqtl_prepare_expression

mkdir /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/02_prep_inputs/snp_annot

python /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/PredictDB-Tutorial/code/split_snp_annot_by_chr.py \
  /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/02_prep_inputs/snp_annot_prep.txt \
  /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/02_prep_inputs/snp_annot/snp_annot

# Final preprocessing step: process the Normalized Expression File
# gunzip /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/goesHyde_mdd_rnaseq_Amygdala.expression.bed.gz
#
# tail -c +2 /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/goesHyde_mdd_rnaseq_Amygdala.expression.bed \
#   > /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/goesHyde_mdd_rnaseq_Amygdala.expression_UNBED.txt

Rscript /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/code/02_prep_inputs/04_process_norm_exp.R

echo "**** Job ends ****"
date
