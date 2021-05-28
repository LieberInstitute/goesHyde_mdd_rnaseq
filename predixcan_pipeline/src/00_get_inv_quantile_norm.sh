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

mkdir processed-data/
mkdir plots/
mkdir code/

ml htslib

module use /jhpce/shared/jhpce/modulefiles/libd

Rscript 00_prepare_gene_expression.R

vcf="../genotype_data/topmed_mdd_602sample_090120_maf005.vcf.gz"
tpm_gct="../processed-data/tpm.gct"
counts_gct="../processed-data/gene_counts.gct"
annotation_gtf="dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf"

tabix --list-chroms ${vcf} > "../processed-data/vcf_chr_list.txt"

# adding normalization script from gtex-pipeline: https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl
wget https://raw.githubusercontent.com/broadinstitute/gtex-pipeline/master/qtl/src/eqtl_prepare_expression.py

eqtl_prepare_expression.py ${tpm_gct} ${counts_gct} ${annotation_gtf} \
    ${sample_participant_lookup} vcf_chr_list.txt ${prefix} \
    --tpm_threshold 0.1 \
    --count_threshold 6 \
    --sample_frac_threshold 0.2 \
    --normalization_method tmm

echo "**** Job ends ****"
date
