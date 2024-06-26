#!/bin/bash
#$ -cwd
#$ -l mem_free=25G,h_vmem=25G,h_fsize=100G
#$ -N generate_ratio_matrix_chrom
#$ -o logs/generate_ratio_matrix_chrom.txt
#$ -e logs/generate_ratio_matrix_chrom.txt

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

module load python/2.7

# Create a ratio file with only autosomes
# echo "
# grep autosomes..."
# zgrep -f ../data/tokeep.txt ../data/clusters/leafcutter_perind.counts.gz > ../data/clusters/leafcutter_perind.counts.chrom

# zip autosome ratio file 
# echo "zip autozome ratio file..."
# gzip ../data/clusters/leafcutter_perind.counts.chrom

# Create qqnormed file for each autosome
echo "build qqnorm file..."
python prepare_phenotype_table.py  ../data/clusters/leafcutter_perind.counts.chrom.gz -p 10

# Compile each autosome into a single file
# head -1 ../data/clusters/leafcutter_perind.counts.chrom.gz.qqnorm_chr1 > ../data/clusters/leafcutter_perind.counts.chrom.gz.qqnorm_all
# cat leafcutter_perind.counts.chrom.gz.qqnorm_chr* | bedtools sort -i - >> leafcutter_perind.counts.chrom.gz.qqnorm_all

echo "**** Job ends ****"
date