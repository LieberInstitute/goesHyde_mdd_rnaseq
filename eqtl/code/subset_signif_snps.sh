#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=12G,h_vmem=12G,h_fsize=200G
#$ -N subset_signif_snps
#$ -o logs/subset_signif_snps.txt
#$ -e logs/subset_signif_snps.txt


module load bcftools



## Locate file and ids

bcftools view --include 'ID=@../data/signif_snps/significant_snps.txt' /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/mdd_bpd/maf01/mdd_bpd_01.vcf.gz > ../data/signif_snps/LIBD_Brain_merged_maf_005_topmed_051120_mdd_signif.vcf

