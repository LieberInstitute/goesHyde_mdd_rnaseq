Deconvolution
========
** Compute and Process **
[`deconvo_data_prep.R`](deconvo_data_prep.R): Run MuSiC analysis, save outputs as .RData
[`deconvo_heatmap.R`](deconvo_heatmap.R)

** Plot helper functions **
[`big_little_boxplot.R`](big_little_boxplot.R): define functions for creating boxplots
divided by scale, makes small proptions visible on plots
[`cell_colors.R`](cell_colors.R): define color scheme for the cell types (including multi for top40 annotation)

** Plot 
[`deconvo_plots.R`](deconvo_plots.R): Create scatterplots comparing datatypes used in MuSiC & boxplots of cell type prop vs. Dx 
[`deconvo_vs_qSV.R`](deconvo_vs_qSV.R): Compare cell type propotions and qSVs for each sample, calculate p-values to detirmine signficance.
[`music_deconvo.R`](music_deconvo.R): Create heatmaps for top40 genes in single cell and bulk data
