#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -N eqtl_prepare_expression
#$ -j y
#$ -o eqtl_prepare_expression.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${TASK_ID}"

cd "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/"

mkdir -p processed-data/01_get_inv_quantile_norm
mkdir plots/
mkdir -p code/01_get_inv_quantile_norm

ml htslib
ml bcftools
ml python/3.6.9

module use /jhpce/shared/jhpce/modulefiles/libd

vcf="../genotype_data/topmed_mdd_602sample_090120_maf005.vcf.gz"
tpm_gct="processed-data/tpm.gct"
counts_gct="processed-data/gene_counts.gct"
annotation_gtf="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf"
prefix="goesHyde_mdd_rnaseq_"

bcftools query -l ${vcf} > "processed-data/01_get_inv_quantile_norm/vcf_samples.txt"

Rscript code/01_get_inv_quantile_norm/01_prepare_gene_expression.R

tabix --list-chroms ${vcf} > "processed-data/01_get_inv_quantile_norm/vcf_chr_list.txt"

# adding normalization script from gtex-pipeline: https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl
wget -P "code/01_get_inv_quantile_norm/" https://raw.githubusercontent.com/broadinstitute/gtex-pipeline/master/qtl/src/eqtl_prepare_expression.py

eqtl_prepare_expression.py ${tpm_gct} ${counts_gct} ${annotation_gtf} \
    ${sample_participant_lookup} "processed-data/vcf_chr_list.txt" ${prefix} \
    --tpm_threshold 0.1 \
    --count_threshold 6 \
    --sample_frac_threshold 0.2 \
    --normalization_method tmm

mv "${prefix}*" "processed-data/01_get_inv_quantile_norm/"

echo "**** Job ends ****"
date
