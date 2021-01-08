#!/bin/bash

## deconvolution directory
MAINDIR="/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/deconvolution"
cd ${MAINDIR}

## Filter sce data 
rm logs/sce_data_prep.txt
qsub sce_data_prep.sh

## Find top marker genes
rm logs/find_markers.txt
qsub find_markers.sh

## Create heatmaps of selected markers
rm logs/deconvo_heatmap.txt
qsub deconvo_heatmap.sh

## Run MuSiC
rm logs/music_deconvo.txt
qsub music_deconvo.sh

## Create plot of deconvolution results
rm logs/deconvo_plots.txt
qsub deconvo_plots.sh

## Compare deconvolution results and qSVs 
rm logs/deconvo_vs_qSV.txt
qsub deconvo_vs_qSV.sh


