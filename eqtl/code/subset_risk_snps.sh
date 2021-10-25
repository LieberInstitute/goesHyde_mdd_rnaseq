#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=12G,h_vmem=12G,h_fsize=200G
#$ -N getVcfRS
#$ -pe local 4
#$ -o logs/getVcfRS$TASK_ID.txt
#$ -e logs/getVcfRS.$TASK_ID.txt
#$ -t 6
#$ -tc 6
#$ -m a





module load bcftools



## Locate file and ids
rs=$(awk 'BEGIN {FS="\t"} {print $1}' ../data/risk_snps/MDD_risk.txt | awk "NR==${SGE_TASK_ID}")



bcftools view -i RS=$rs /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/mdd_bpd/maf01/mdd_bpd_01.vcf.gz > ../data/risk_snps/MDD_temp/mdd_bpd_01_$rs.vcf
