deconvolution
========

[`deconvolution.sh`](deconvolution.sh) Runs full deconvolution analysis

### Prep sce data & Find Markers
1. [`sce_data_prep.R`](sce_data_prep.R): Filter cell types, filter for common genes
with bulk data sets and for genes with all 0 expression , add bp lenghts (for RPKM, not used for now) 
2. [`find_markers.R`](find_markers.R): Run [`get_mean_ratio()`](get_mean_ratio.R)
 and [`findMarkers_1vAll()`](findMarkers_1vAll.R) on the sce data to get produce 
 stats for all genes. Save marker stats as Rdata and top 5 genes as csv. Produce 
 expression and mean_ratio vs. std.logFC plots.
3. [`deconvo_heatmap.R`](deconvo_heatmap.R): Create heatmaps for selected marker genes in 
single cell, bulk data, and external Velmechev sce data. This is to confirm that 
"clean" markers have been selected

### Run MuSiC & Plot results
4. [`music_deconvo.R`](music_deconvo.R): Run MuSiC analysis on all combonations of
cellTypes (specific & broad) and regions, save outputs and  as .RData
5. [`deconvo_plots.R`](deconvo_plots.R): Create boxplots of the MuSiC output
6. [`deconvo_vs_qSV.R`](deconvo_vs_qSV.R): Compare cell type propotions and qSVs 
for each sample, calculate p-values to detirmine signficance.

### R functions 
[`get_mean_ratio.R`](get_mean_ratio.R): Filter genes for median expresssion in 
target cell type != 0, calculate ratio of mean expession of target genes/highest
mean expession from non-target cell type.

[`findMarkers_1vAll.R`](findMarkers_1vAll.R) Apply findMarker in 1 vs. All fashion.

[`big_little_boxplot.R`](big_little_boxplot.R): define functions for creating 
boxplots divided by scale, makes small proptions visible on plots

[`cell_colors.R`](cell_colors.R): define color scheme for the cell types 
(including multi annotation)


