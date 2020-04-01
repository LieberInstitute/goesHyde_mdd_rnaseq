#!/bin/bash
#$ -cwd
#$ -pe local 8
#$ -l mem_free=28G,h_vmem=30G,h_fsize=200G
#$ -N step7-Rcounts-goesHyde.MDD
#$ -o ./logs/Rcounts-goesHyde.txt
#$ -e ./logs/Rcounts-goesHyde.txt
#$ -hold_jid pipeline_setup,step4-featCounts-goesHyde.MDD,step6-txQuant-goesHyde.MDD
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

## Force R 3.3.x in JHPCE (to avoid some issues with conda_R)
module unload conda_R
module load R/3.3.x

Rscript /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/.step7-create_count_objects-human.R -o hg38 -m /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data -e goesHyde -p MDD -l TRUE -c TRUE -t 8 -s reverse

echo "**** Job ends ****"
date
