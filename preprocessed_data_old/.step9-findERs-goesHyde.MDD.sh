#!/bin/bash
#$ -cwd
#$ -l mem_free=25G,h_vmem=25G,h_fsize=200G
#$ -N step9-findERs-goesHyde.MDD
#$ -o ./logs/findERs-goesHyde.txt
#$ -e ./logs/findERs-goesHyde.txt
#$ -hold_jid pipeline_setup,step5b-meanCoverage-goesHyde.MDD
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

for meanFile in Coverage/mean*.bw
do
    echo "************************************"
    date
    echo "Initializing script for ${meanFile}"
    echo "************************************"
    Rscript /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/.step9-find_expressed_regions.R -m ${meanFile} -o /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/ERs -i /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.gencode
done

echo "**** Job ends ****"
date
