#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -N eqtl_prepare_expression
#$ -j y
#$ -o /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/logs/eqtl_prepare_expression.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${TASK_ID}"

cd "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/"

mkdir -p processed-data/01_get_inv_quantile_norm
mkdir plots/

ml htslib
ml bcftools
ml plink/1.90b6.6
ml conda_R/4.1
ml samtools/1.10

module use /jhpce/shared/jhpce/modulefiles/libd
# TODO vcf might not be filtered enough - Hardy-Weinberg, 0.1%, too many SNPs, may need to filter eventually
# Make sure if they recommend using SNPs not in LD - may prune out SNPs in LD

# Sample IDs in the VCF, GCTs and sample lookup MUST BE IN SAME ORDER
Amyg_tpm_gct="/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/Amygdala_tpm.gct"
Amyg_counts_gct="/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/Amygdala_gene_counts.gct"

sACC_tpm_gct="/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/sACC_tpm.gct"
sACC_counts_gct="/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/sACC_gene_counts.gct"

annotation_gtf="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf"
vcf_chr_list="/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/vcf_chr_list.txt"

Amyg_sample_participant_lookup="/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/Amygdala_samp_part_lookup.txt"

Rscript code/01_get_inv_quantile_norm/01_prepare_gene_expression.R -r "Amygdala"

Amyg_vcf="/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/topmed_mdd_602sample_090120_maf005_Amygdala_sorted.vcf.gz"
# sACC_vcf="processed-data/01_get_inv_quantile_norm/topmed_mdd_602sample_090120_maf005_sACC_sorted.vcf.gz"

prefix_Amygdala="goesHyde_mdd_rnaseq_Amygdala"

bgzip -c "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/topmed_mdd_602sample_090120_maf005_Amygdala_sorted.vcf" > ${Amyg_vcf}
# bgzip -c "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/topmed_mdd_602sample_090120_maf005_sACC_sorted.vcf" > ${sACC_vcf}

tabix -p vcf ${Amyg_vcf}
tabix -p vcf ${sACC_vcf}

# It's just 1-23
tabix --list-chroms ${Amyg_vcf} | awk '$1="chr"$1' > "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/vcf_chr_list.txt"

# if [ ! -f "code/01_get_inv_quantile_norm/eqtl_prepare_expression.py" ]; then
#   # adding normalization script from gtex-pipeline: https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl
#   wget -P "code/01_get_inv_quantile_norm/" https://raw.githubusercontent.com/broadinstitute/gtex-pipeline/master/qtl/src/eqtl_prepare_expression.py
# else
#   echo "eqtl_prepare_expression.py already exists."
# fi

# wget -O KJ_eqtl_prepare_expression.py https://raw.githubusercontent.com/KrotosBenjamin/gtex-pipeline/master/qtl/src/eqtl_prepare_expression.py
# git clone git@github.com:francois-a/rnaseq-utils.git

# git submodule add https://github.com/francois-a/rnaseq-utils.git

# conda create -n eqtl_prepare_expression python=3.5
# pip install pandas numpy scipy argparse qtl --user
conda activate eqtl_prepare_expression
export PYTHONPATH=$PYTHONPATH:/$PWD/rnaseq-utils

# bgzip needs to be mounted as a module through htslib for this step. If you're
# encountering some bgzip related error, it's probably for that reason.
# python -m pdb
# b 118

# https://groups.google.com/g/predixcanmetaxcan/c/iC-aXrMUdT8/m/Z71Dgl1BBwAJ
# "PredictDB-tutorial (https://github.com/hakyimlab/PredictDB-Tutorial) might
# help you once you have normalized expression (via any normalization of your choice)."
./code/01_get_inv_quantile_norm/KJ_eqtl_prepare_expression.py ${Amyg_tpm_gct} ${Amyg_counts_gct} ${annotation_gtf} \
    ${Amyg_sample_participant_lookup} ${vcf_chr_list} ${prefix_Amygdala} \
    --tpm_threshold 0.1 \
    --count_threshold 6 \
    --sample_frac_threshold 0.2 \
    --normalization_method tmm \
    --output_dir "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/"

echo "**** Job ends ****"
date
