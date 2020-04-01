#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -N step0-ercc-goesHyde.MDD
#$ -pe local 5
#$ -o logs/ercc-goesHyde.$TASK_ID.txt
#$ -e logs/ercc-goesHyde.$TASK_ID.txt
#$ -t 1-634
#$ -tc 20
#$ -hold_jid pipeline_setup,step00-merge-goesHyde.MDD
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

## Locate file and ids
FILE1=$(awk 'BEGIN {FS="\t"} {print $1}' samples.manifest | awk "NR==${SGE_TASK_ID}")
if [ TRUE == "TRUE" ] 
then
    FILE2=$(awk 'BEGIN {FS="\t"} {print $3}' samples.manifest | awk "NR==${SGE_TASK_ID}")
fi
ID=$(cat samples.manifest | awk '{print $NF}' | awk "NR==${SGE_TASK_ID}")
mkdir -p /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/Ercc/${ID}

if [ TRUE == "TRUE" ]
then 
    /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/kallisto quant     -i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ERCC/ERCC92.idx     -o /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/Ercc/${ID} -t 5 --rf-stranded     ${FILE1} ${FILE2}
else
    /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/kallisto quant     -i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/ERCC/ERCC92.idx     -o /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/Ercc/${ID} -t 5 --single --rf-stranded ${FILE1}
fi

echo "**** Job ends ****"
date
