#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=2G,h_vmem=2G,h_fsize=100G
#$ -N step8-old
#$ -o logs/callVariants-goesHyde.$TASK_ID.txt
#$ -e logs/callVariants-goesHyde.$TASK_ID.txt
#$ -t 540-634
#$ -tc 100
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"
echo "Sample id: $(cat samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")"
echo "****"

mkdir -p /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data_old/Genotypes/

ID=$(cat samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")
BAM=/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data_old/HISAT2_out/${ID}_accepted_hits.sorted.bam

SNPTMP=/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data_old/Genotypes/${ID}_tmp.vcf
SNPOUTGZ=/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data_old/Genotypes/${ID}.vcf.gz

module load bcftools
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/samtools-1.2/samtools mpileup -l /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Genotyping/common_missense_SNVs_hg38.bed -AB -q0 -Q13 -d1000000 -uf /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/GRCh38.primary_assembly.genome.fa ${BAM} -o ${SNPTMP}
bcftools call -mv -Oz ${SNPTMP} > ${SNPOUTGZ}

module load htslib
tabix -p vcf ${SNPOUTGZ}

rm ${SNPTMP}

echo "**** Job ends ****"
date
