#!/bin/bash

## To save typing space later
MAINDIR="/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq"

## Clean the data
cd ${MAINDIR}/data
rm logs/clean_data.txt ## wipe out the previous log since we only care about the last version
qsub clean_data.sh

## Filter the expression features
cd ${MAINDIR}/exprs_cutoff
rm logs/get_expression_cutoff.txt
qsub get_expression_cutoff.sh

## Get degredation regions
cd ${MAINDIR}
rm logs/get_degredation_regions.sh
qsub get_degredation_regions.sh

## Merge MDDseq and BiPseq data
cd ${MAINDIR}/data
rm logs/merge_MDDseq_BiPseq_degradation.txt
qsub merge_MDDseq_BiPseq_degradation.sh

## Run qSV DE analysis
cd ${MAINDIR}/differential_expression
rm logs/qSV_model_DE_analysis.txt
qsub qSV_model_DE_analysis.sh

## Run wgna scripts
cd ${MAINDIR}/wgna
rm logs/run_wgcna_combined.txt
qsub run_wgcna_combined.sh

