#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=12G,h_vmem=12G,h_fsize=200G
#$ -N subset_signif_snps
#$ -o logs/subset_signif_snps_MDD.txt
#$ -e logs/subset_signif_snps_MDD.txt


module load bcftools
module load plink


## Locate file and ids

bcftools view --include 'ID=@../data/risk_snps/Variants_within_MDD_loci_all_052722.txt' /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/mdd_bpd/maf01/mdd_bpd_01.vcf.gz > ../data/risk_snps/LIBD_maf01_gwas_MDD.vcf.gz
plink2 --vcf ../data/risk_snps/LIBD_maf01_gwas_MDD.vcf.gz --make-bed --out ../data/risk_snps/LIBD_maf01_gwas_MDD

