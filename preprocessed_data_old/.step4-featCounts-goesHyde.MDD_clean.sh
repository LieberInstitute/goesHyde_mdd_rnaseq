#!/bin/bash
#$ -cwd
#$ -N step4-featCounts-goesHyde.MDD_clean
#$ -o ./logs/featCounts-goesHyde_clean.txt
#$ -e ./logs/featCounts-goesHyde_clean.txt
#$ -hold_jid pipeline_setup,step4-featCounts-goesHyde.MDD
#$ -m a
echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "****"

## Delete temporary files after they have been used
rm -rf /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/Counts/junction/tmpdir

echo "**** Job ends ****"
date
