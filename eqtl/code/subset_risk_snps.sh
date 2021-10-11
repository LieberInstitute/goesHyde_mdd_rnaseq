#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=12G,h_vmem=12G,h_fsize=200G
#$ -N getVcfRS
#$ -pe local 4
#$ -o logs/getVcfRS$TASK_IDo.txt
#$ -e logs/getVcfRS.$TASK_IDe.txt
#$ -t 6
#$ -tc 6
#$ -m a





module load bcftools



## Locate file and ids
rs=$(awk 'BEGIN {FS="\t"} {print $1}' ../data/risk_snps/MDD_risk_test.txt | awk "NR==${SGE_TASK_ID}")



bcftools view -i RS=$rs /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/mdd_bpd/maf01/mdd_bpd_01.vcf.gz > ../data/risk_snps/LIBD_Brain_merged_maf_005_topmed_051120_mdd_$rs.vcf
