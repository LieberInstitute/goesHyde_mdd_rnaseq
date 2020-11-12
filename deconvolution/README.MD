deconvolution
========

## Compute and Process

[`music_deconvo.R`](music_deconvo.R): Run MuSiC analysis, save outputs as .RData
[`deconvo_data_prep.R`](deconvo_data_prep.R): Make long versions of prop data, 
add Dx and age from res_gene object

## Plot helper functions

[`big_little_boxplot.R`](big_little_boxplot.R): define functions for creating 
boxplots divided by scale, makes small proptions visible on plots
[`cell_colors.R`](cell_colors.R): define color scheme for the cell types 
(including multi for top40 annotation)

## Plot

[`deconvo_plots.R`](deconvo_plots.R): Create scatterplots comparing datatypes 
used in MuSiC & boxplots of cell type prop vs. Dx 
[`deconvo_vs_qSV.R`](deconvo_vs_qSV.R): Compare cell type propotions and qSVs 
for each sample, calculate p-values to detirmine signficance.
[`deconvo_heatmap.R`](deconvo_heatmap.R): Create heatmaps for top40 genes in 
single cell and bulk data
