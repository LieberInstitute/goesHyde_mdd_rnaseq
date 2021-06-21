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

ml htslib
ml bcftools
ml conda_R/4.1

module use /jhpce/shared/jhpce/modulefiles/libd
# TODO vcf might not be filtered enough - Hardy-Weinberg, 0.1%, too many SNPs, may need to filter eventually
# Make sure if they recommend using SNPs not in LD - may prune out SNPs in LD

# Sample IDs in the VCF, GCTs and sample lookup MUST BE IN SAME ORDER
Amyg_tpm_gct="processed-data/01_get_inv_quantile_norm/Amygdala_tpm.gct"
Amyg_counts_gct="processed-data/01_get_inv_quantile_norm/Amygdala_gene_counts.gct"

sACC_tpm_gct="processed-data/01_get_inv_quantile_norm/sACC_tpm.gct"
sACC_counts_gct="processed-data/01_get_inv_quantile_norm/sACC_gene_counts.gct"

annotation_gtf="/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf"
prefix="goesHyde_mdd_rnaseq_"
vcf_chr_list="processed-data/01_get_inv_quantile_norm/vcf_chr_list.txt"

Amyg_sample_participant_lookup="processed-data/01_get_inv_quantile_norm/Amygdala_samp_part_lookup.txt"

Rscript code/01_get_inv_quantile_norm/01_prepare_gene_expression.R

Amyg_vcf="processed-data/01_get_inv_quantile_norm/topmed_mdd_602sample_090120_maf005_Amygdala_sorted.vcf.gz"
sACC_vcf="processed-data/01_get_inv_quantile_norm/topmed_mdd_602sample_090120_maf005_sACC_sorted.vcf.gz"

bgzip -c "processed-data/01_get_inv_quantile_norm/topmed_mdd_602sample_090120_maf005_Amygdala_sorted.vcf" > ${Amyg_vcf}
bgzip -c "processed-data/01_get_inv_quantile_norm/topmed_mdd_602sample_090120_maf005_sACC_sorted.vcf" > ${sACC_vcf}

tabix -p vcf ${Amyg_vcf}
tabix -p vcf ${sACC_vcf}

# It's just 1-23
tabix --list-chroms ${Amyg_vcf} > "processed-data/01_get_inv_quantile_norm/vcf_chr_list.txt"

if [ ! -f "code/01_get_inv_quantile_norm/eqtl_prepare_expression.py" ]; then
  # adding normalization script from gtex-pipeline: https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl
  wget -P "code/01_get_inv_quantile_norm/" https://raw.githubusercontent.com/broadinstitute/gtex-pipeline/master/qtl/src/eqtl_prepare_expression.py
else
  echo "eqtl_prepare_expression.py already exists."
fi

# conda create -n eqtl_prepare_expression python=3.5
# pip install pandas numpy scipy argparse qtl --user
conda activate eqtl_prepare_expression

./code/01_get_inv_quantile_norm/eqtl_prepare_expression.py ${Amyg_tpm_gct} ${Amyg_counts_gct} ${annotation_gtf} \
    ${Amyg_sample_participant_lookup} ${vcf_chr_list} "goesHyde_mdd_rnaseq_Amygdala" \
    --tpm_threshold 0.1 \
    --count_threshold 6 \
    --sample_frac_threshold 0.2 \
    --normalization_method tmm

# 16:13 predixcan_pipeline $ ./code/01_get_inv_quantile_norm/eqtl_prepare_expression.py ${tpm_gct} ${counts_gct} ${annotation_gtf} \
# >     ${sample_participant_lookup} ${vcf_chr_list} ${prefix} \
# >     --tpm_threshold 0.1 \
# >     --count_threshold 6 \
# >     --sample_frac_threshold 0.2 \
# >     --normalization_method tmm
# Loading expression data
# sys:1: DtypeWarning: Columns (15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326) have mixed types. Specify dtype option on import or set low_memory=False.
# Traceback (most recent call last):
#   File "./code/01_get_inv_quantile_norm/eqtl_prepare_expression.py", line 97, in <module>
#     raise ValueError('Sample IDs in expression files and participant lookup table must match.')
# ValueError: Sample IDs in expression files and participant lookup table must match.

mv "${prefix}*" "processed-data/01_get_inv_quantile_norm/"

echo "**** Job ends ****"
date

# mv code/01_get_inv_quantile_norm eqtl_prepare_expression.log processed-data/01_get_inv_quantile_norm
