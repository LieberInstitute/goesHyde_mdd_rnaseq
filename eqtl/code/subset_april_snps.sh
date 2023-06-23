#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=12G,h_vmem=12G,h_fsize=200G
#$ -N subset_signif_snps
#$ -o logs/subset_signif_snps.txt
#$ -e logs/subset_signif_snps.txt


module load bcftools



## Locate file and ids

bcftools view --include 'ID=@../data/plot_data/april_snps.txt' /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/mdd_bpd/maf01/mdd_bpd_01.vcf.gz > ../data/plot_data/LIBD_Brain_merged_maf_005_topmed_april_snps.vcf.gz

