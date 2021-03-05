#!/bin/bash

## These mkdir steps + ln -s + "mkdir -p logs/amygdala_gene" were typically done
## outside the loop at
## https://github.com/LieberInstitute/twas/blob/master/bsp2/compute_weights_indv.sh

## To avoid having to change the file permissions later
## From https://twitter.com/fellgernon/status/1258455434073124865?s=20
## and https://www.cyberciti.biz/tips/understanding-linux-unix-umask-value-usage.html
umask u=rwx,g=rwx,o= ## equivalent to umask 007

## Required order for running this code:

# For the logs
rm -r trash/sACC_gene
mv logs/sACC_gene trash/
mkdir -p logs/sACC_gene

## For output files
mkdir -p sACC_gene/tmp_files
mkdir -p sACC_gene/out_files

# For GEMMA
ln -s /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/twas/sACC_gene/ sACC_gene/output

## For running the main script
qsub compute_weights_indv_sACC_full_gene.sh