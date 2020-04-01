#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -N step3-hisat2
#$ -pe local 8
#$ -o logs/hisat2-test.txt
#$ -e logs/hisat2-test.txt
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"
echo "****"

# Directories
mkdir -p /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/HISAT2_out/test

## Locate file and ids
FILE1="/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/FASTQ_merged/R14179.fastq.gz"
FILE2="/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/FASTQ_merged/R14179_read2.fastq.gz"
ID="R14179"

echo "HISAT2 alignment run on trimmed paired-end reads"
FP=/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/trimmed_fq/${ID}_trimmed_forward_paired.fastq
FU=/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/trimmed_fq/${ID}_trimmed_forward_unpaired.fastq
RP=/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/trimmed_fq/${ID}_trimmed_reverse_paired.fastq
RU=/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/trimmed_fq/${ID}_trimmed_reverse_unpaired.fastq

/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/hisat2-2.0.4/hisat2 -p 8 	-x /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/hisat2_GRCh38primary -1 $FP -2 $RP -U ${FU},${RU} 	-S /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/HISAT2_out/test/${ID}_hisat_out1.sam --rna-strandness RF --phred33      	2>/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/HISAT2_out/test/${ID}_trimmed_summary.txt

## Untrimmed, pair-end
echo "HISAT2 alignment run on original untrimmed paired-end reads"
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/hisat2-2.0.4/hisat2 -p 8 	-x /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/hisat2_GRCh38primary -1 ${FILE1} -2 ${FILE2} 	-S /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/HISAT2_out/test/${ID}_hisat_out2.sam --rna-strandness RF --phred33      	2>/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/HISAT2_out/test/${ID}_original_summary.txt



## Locate file and ids
FILE1="/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/FASTQ_merged/R17609.fastq.gz"
FILE2="/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/FASTQ_merged/R17609_read2.fastq.gz"
ID="R17609"

echo "HISAT2 alignment run on trimmed paired-end reads"
FP=/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/trimmed_fq/${ID}_trimmed_forward_paired.fastq
FU=/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/trimmed_fq/${ID}_trimmed_forward_unpaired.fastq
RP=/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/trimmed_fq/${ID}_trimmed_reverse_paired.fastq
RU=/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/trimmed_fq/${ID}_trimmed_reverse_unpaired.fastq

/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/hisat2-2.0.4/hisat2 -p 8 	-x /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/hisat2_GRCh38primary -1 $FP -2 $RP -U ${FU},${RU} 	-S /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/HISAT2_out/test/${ID}_hisat_out1.sam --rna-strandness RF --phred33      	2>/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/HISAT2_out/test/${ID}_trimmed_summary.txt

## Untrimmed, pair-end
echo "HISAT2 alignment run on original untrimmed paired-end reads"
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/hisat2-2.0.4/hisat2 -p 8 	-x /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/hisat2_GRCh38primary -1 ${FILE1} -2 ${FILE2} 	-S /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/HISAT2_out/test/${ID}_hisat_out2.sam --rna-strandness RF --phred33      	2>/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/HISAT2_out/test/${ID}_original_summary.txt

echo "**** Job ends ****"
date
