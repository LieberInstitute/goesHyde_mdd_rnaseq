#!/bin/bash
#$ -cwd
#$ -l mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -N bed_bgzip
#$ -o logs/02_bed_bgzip.txt
#$ -e logs/02_bed_bgzip.txt
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
module load htslib

## List current modules for reproducibility
module list

## Edit with your job command
bgzip /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_input/expression_bed/gene_Amygdala.bed && tabix -p bed /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_input/expression_bed/gene_Amygdala.bed.gz
bgzip /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_input/expression_bed/gene_sACC.bed && tabix -p bed /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_input/expression_bed/gene_sACC.bed.gz
bgzip /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_input/expression_bed/exon_Amygdala.bed && tabix -p bed /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_input/expression_bed/exon_Amygdala.bed.gz
bgzip /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_input/expression_bed/exon_sACC.bed && tabix -p bed /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_input/expression_bed/exon_sACC.bed.gz
bgzip /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_input/expression_bed/jxn_Amygdala.bed && tabix -p bed /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_input/expression_bed/jxn_Amygdala.bed.gz
bgzip /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_input/expression_bed/jxn_sACC.bed && tabix -p bed /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_input/expression_bed/jxn_sACC.bed.gz
bgzip /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_input/expression_bed/tx_Amygdala.bed && tabix -p bed /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_input/expression_bed/tx_Amygdala.bed.gz
bgzip /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_input/expression_bed/tx_sACC.bed && tabix -p bed /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_input/expression_bed/tx_sACC.bed.gz

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/