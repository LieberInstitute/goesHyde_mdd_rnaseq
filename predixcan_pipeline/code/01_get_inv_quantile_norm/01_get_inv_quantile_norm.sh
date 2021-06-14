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
ml conda_R/4.1

module use /jhpce/shared/jhpce/modulefiles/libd
# TODO vcf might not be filtered enough - Hardy-Weinberg, 0.1%, too many SNPs, may need to filter eventually
# Make sure if they recommend using SNPs not in LD - may prune out SNPs in LD

# Sample IDs in the VCF, GCTs and sample lookup MUST BE IN SAME ORDER
tpm_gct="processed-data/01_get_inv_quantile_norm/tpm.gct"
counts_gct="processed-data/01_get_inv_quantile_norm/gene_counts.gct"
annotation_gtf="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf"
prefix="goesHyde_mdd_rnaseq_"
vcf_chr_list="processed-data/01_get_inv_quantile_norm/vcf_chr_list.txt"
sample_participant_lookup="processed-data/01_get_inv_quantile_norm/samp_part_lookup.txt"

Rscript code/01_get_inv_quantile_norm/01_prepare_gene_expression.R

vcf="processed-data/01_get_inv_quantile_norm/topmed_mdd_602sample_090120_maf005_MDD_sorted.vcf" # get more info about this - ask Josh?

bgzip -c ${vcf} > "processed-data/01_get_inv_quantile_norm/topmed_mdd_602sample_090120_maf005_MDD_sorted.vcf.gz"
vcf="processed-data/01_get_inv_quantile_norm/topmed_mdd_602sample_090120_maf005_MDD_sorted.vcf.gz"
tabix -p vcf ${vcf}
tabix --list-chroms ${vcf} > "processed-data/01_get_inv_quantile_norm/vcf_chr_list.txt"

if [ ! -f "code/01_get_inv_quantile_norm/eqtl_prepare_expression.py" ]; then
  # adding normalization script from gtex-pipeline: https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl
  wget -P "code/01_get_inv_quantile_norm/" https://raw.githubusercontent.com/broadinstitute/gtex-pipeline/master/qtl/src/eqtl_prepare_expression.py
else
  echo "eqtl_prepare_expression.py already exists."
fi

# conda create -n eqtl_prepare_expression python=3.5
# pip install pandas numpy scipy argparse qtl --user
conda activate eqtl_prepare_expression

./code/01_get_inv_quantile_norm/eqtl_prepare_expression.py ${tpm_gct} ${counts_gct} ${annotation_gtf} \
    ${sample_participant_lookup} ${vcf_chr_list} ${prefix} \
    --tpm_threshold 0.1 \
    --count_threshold 6 \
    --sample_frac_threshold 0.2 \
    --normalization_method tmm

mv "${prefix}*" "processed-data/01_get_inv_quantile_norm/"

echo "**** Job ends ****"
date

# mv code/01_get_inv_quantile_norm eqtl_prepare_expression.log processed-data/01_get_inv_quantile_norm
