
library(SummarizedExperiment)
library(SingleCellExperiment)
library(jaffelab)
library(pheatmap)
library(DeconvoBuddies)
library(tidyverse)
library(RColorBrewer)
library(rlang)
library(recount)
library(here)
library(sessioninfo)

# Load colors
source(here("main_colors.R"))
load(here("deconvolution","data","cell_colors.Rdata"), verbose = TRUE)

## Load rse_gene data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

#### Marker genes ####
load(here("deconvolution","data","marker_genes.Rdata"), verbose = TRUE)
load(here("deconvolution","data","marker_stats.Rdata"), verbose = TRUE)

#### sce data ####
load(here("deconvolution","data","sce_filtered.Rdata"), verbose = TRUE)
sce_filter <- sce_pan[marker_genes,]
dim(sce_filter)
# [1]   150 34005

pb <- pseudobulk(sce_filter,  cell_group_cols = c("cellType.Broad","donor","region"), add_symbol = TRUE)
dim(pb)
# [1] 150  70

#### Prep for heatmaps ####
## define annotation tables

pb_anno <- data.frame(names = colnames(pb), names2 = colnames(pb)) %>%
  column_to_rownames("names") %>%
  separate(names2, into = c("cellType.Broad","donor","Region"))

## create gene anno obj
marker_anno <- marker_stats %>%
  filter(rank_ratio <= 25) %>%
  select(gene, Symbol, cellType.target, rank_ratio) 

marker_anno_gene <- marker_anno %>%
  select(-Symbol) %>%
  column_to_rownames(var = "gene")
                     
marker_anno_symb <- marker_anno %>%
  select(-gene) %>%
  column_to_rownames(var = "Symbol")

## define donor colrs
donors <- unique(pb_anno$donor)

donor_colors <- brewer.pal(n= length(donors), name = "Greys")
names(donor_colors) <- donors
cell_colors <- cell_colors[levels(sce_filter$cellType.Broad)]
anno_colors <- list(cellType.target = cell_colors,
                    cellType.Broad = cell_colors,
                    donor = donor_colors,
                    PrimaryDx = mdd_Dx_colors,
                    Experiment = mdd_dataset_colors,
                    BrainRegion = mdd_BrainRegion_colors[c("Amygdala","sACC")])

#### Create Heatmaps SCE ####
png(here("deconvolution","plots","heatmap_markers-sce.png"), height = 750, width = 550)
pheatmap(pb,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_row = marker_anno_symb,
         annotation_col = pb_anno,
         annotation_colors = anno_colors,
         cluster_rows = TRUE,
         main = "Pan-brain data:log2(pseudobulk counts) - Top 25 Mean Ratio"
)
dev.off()


png(here("deconvolution","plots","heatmap_markers-sce_long.png"), height = 2000, width = 550)
pheatmap(pb,
         show_colnames = FALSE,
         show_rownames = TRUE,
         annotation_row = marker_anno_symb,
         annotation_col = pb_anno,
         annotation_colors = anno_colors,
         cluster_rows = TRUE,
         main = "Pan-brain data:log2(pseudobulk counts) - Top 25 Mean Ratio"
)
dev.off()

png(here("deconvolution","plots","heatmap_markers-sce_unsorted.png"), height = 750, width = 550)
pheatmap(pb,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_row = marker_anno_symb,
         annotation_col = pb_anno,
         annotation_colors = anno_colors,
         cluster_rows = FALSE,
         main = "Pan-brain data:log2(pseudobulk counts) - Top 25 Mean Ratio"
)
dev.off()

#### Bulk Heatmaps ####
bulk_filtered <- rse_gene[marker_genes, ]
dim(bulk_filtered)
bulk_RPKM <- log(recount::getRPKM(bulk_filtered, "Length")+1)
bulk_anno <- as.data.frame(colData(rse_gene)[,c("PrimaryDx",'Experiment',"BrainRegion")])

png(paste0("plots/heatmap_markers-bulk.png"), height = 750, width = 800)
pheatmap(bulk_RPKM,
         show_colnames = FALSE,
         show_rownames = FALSE,
         annotation_row = marker_anno_gene,
         annotation_col = bulk_anno,
         annotation_colors = anno_colors,
         main = "Bulk log(RPKM +1)")
dev.off()

png(paste0("plots/heatmap_markers-bulk_long.png"), height = 2000, width = 800)
bulk_heatmap <- pheatmap(bulk_RPKM,
         show_colnames = FALSE,
         show_rownames = TRUE,
         annotation_row = marker_anno_symb,
         annotation_col = bulk_anno,
         annotation_colors = anno_colors,
         main = "Bulk log(RPKM +1)")
dev.off()

problem_genes <- bulk_heatmap$tree_row$labels[bulk_heatmap$tree_row$order[1:2]]
save(problem_genes, file = here("deconvolution","data","problem_genes.Rdata"))
marker_stats %>% filter(gene %in% problem_genes)

names(bulk_heatmap)



# sgejobs::job_single('deconvo_heatmap.R', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript deconvo_heatmap.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# [1] "2021-01-08 18:03:09 EST"
# user  system elapsed 
# 243.500  17.748 263.964 
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value                                      
# version  R version 4.0.3 Patched (2020-11-29 r79529)
# os       CentOS Linux 7 (Core)                      
# system   x86_64, linux-gnu                          
# ui       X11                                        
# language (EN)                                       
# collate  en_US.UTF-8                                
# ctype    en_US.UTF-8                                
# tz       US/Eastern                                 
# date     2021-01-08                                 
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date       lib source                                   
# AnnotationDbi          1.52.0   2020-10-27 [2] Bioconductor                             
# askpass                1.1      2019-01-13 [2] CRAN (R 4.0.3)                           
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.0.3)                           
# backports              1.2.0    2020-11-02 [1] CRAN (R 4.0.3)                           
# base64enc              0.1-3    2015-07-28 [2] CRAN (R 4.0.3)                           
# beachmat               2.6.2    2020-11-24 [1] Bioconductor                             
# Biobase              * 2.50.0   2020-10-27 [2] Bioconductor                             
# BiocFileCache          1.14.0   2020-10-27 [2] Bioconductor                             
# BiocGenerics         * 0.36.0   2020-10-27 [2] Bioconductor                             
# BiocParallel           1.24.1   2020-11-06 [2] Bioconductor                             
# biomaRt                2.46.0   2020-10-27 [2] Bioconductor                             
# Biostrings             2.58.0   2020-10-27 [2] Bioconductor                             
# bit                    4.0.4    2020-08-04 [2] CRAN (R 4.0.3)                           
# bit64                  4.0.5    2020-08-30 [2] CRAN (R 4.0.3)                           
# bitops                 1.0-6    2013-08-17 [2] CRAN (R 4.0.3)                           
# blob                   1.2.1    2020-01-20 [2] CRAN (R 4.0.3)                           
# broom                  0.7.3    2020-12-16 [2] CRAN (R 4.0.3)                           
# BSgenome               1.58.0   2020-10-27 [2] Bioconductor                             
# bumphunter             1.32.0   2020-10-27 [2] Bioconductor                             
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.0.3)                           
# checkmate              2.0.0    2020-02-06 [2] CRAN (R 4.0.3)                           
# cli                    2.2.0    2020-11-20 [1] CRAN (R 4.0.3)                           
# cluster                2.1.0    2019-06-19 [3] CRAN (R 4.0.3)                           
# codetools              0.2-18   2020-11-04 [3] CRAN (R 4.0.3)                           
# colorspace             2.0-0    2020-11-11 [2] CRAN (R 4.0.3)                           
# crayon                 1.3.4    2017-09-16 [2] CRAN (R 4.0.3)                           
# curl                   4.3      2019-12-02 [2] CRAN (R 4.0.3)                           
# data.table             1.13.6   2020-12-30 [2] CRAN (R 4.0.3)                           
# DBI                    1.1.0    2019-12-15 [2] CRAN (R 4.0.3)                           
# dbplyr                 2.0.0    2020-11-03 [2] CRAN (R 4.0.3)                           
# DelayedArray           0.16.0   2020-10-27 [2] Bioconductor                             
# DelayedMatrixStats     1.12.1   2020-11-24 [1] Bioconductor                             
# derfinder              1.24.2   2020-12-18 [2] Bioconductor                             
# derfinderHelper        1.24.1   2020-12-18 [2] Bioconductor                             
# digest                 0.6.27   2020-10-24 [1] CRAN (R 4.0.3)                           
# doRNG                  1.8.2    2020-01-27 [2] CRAN (R 4.0.3)                           
# downloader             0.4      2015-07-09 [2] CRAN (R 4.0.3)                           
# dplyr                * 1.0.2    2020-08-18 [1] CRAN (R 4.0.3)                           
# ellipsis               0.3.1    2020-05-15 [2] CRAN (R 4.0.3)                           
# fansi                  0.4.1    2020-01-08 [2] CRAN (R 4.0.3)                           
# farver                 2.0.3    2020-01-16 [2] CRAN (R 4.0.3)                           
# forcats              * 0.5.0    2020-03-01 [2] CRAN (R 4.0.3)                           
# foreach                1.5.1    2020-10-15 [2] CRAN (R 4.0.3)                           
# foreign                0.8-81   2020-12-22 [3] CRAN (R 4.0.3)                           
# Formula                1.2-4    2020-10-16 [2] CRAN (R 4.0.3)                           
# fs                     1.5.0    2020-07-31 [1] CRAN (R 4.0.3)                           
# generics               0.1.0    2020-10-31 [2] CRAN (R 4.0.3)                           
# GenomeInfoDb         * 1.26.2   2020-12-08 [2] Bioconductor                             
# GenomeInfoDbData       1.2.4    2020-11-30 [2] Bioconductor                             
# GenomicAlignments      1.26.0   2020-10-27 [2] Bioconductor                             
# GenomicFeatures        1.42.1   2020-11-12 [2] Bioconductor                             
# GenomicFiles           1.26.0   2020-10-27 [2] Bioconductor                             
# GenomicRanges        * 1.42.0   2020-10-27 [2] Bioconductor                             
# GEOquery               2.58.0   2020-10-27 [2] Bioconductor                             
# ggplot2              * 3.3.3    2020-12-30 [2] CRAN (R 4.0.3)                           
# glue                   1.4.2    2020-08-27 [1] CRAN (R 4.0.3)                           
# googledrive            1.0.1    2020-05-05 [1] CRAN (R 4.0.3)                           
# gridExtra              2.3      2017-09-09 [2] CRAN (R 4.0.3)                           
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.0.3)                           
# haven                  2.3.1    2020-06-01 [1] CRAN (R 4.0.3)                           
# here                 * 1.0.0    2020-11-15 [1] CRAN (R 4.0.3)                           
# Hmisc                  4.4-2    2020-11-29 [2] CRAN (R 4.0.3)                           
# hms                    0.5.3    2020-01-08 [2] CRAN (R 4.0.3)                           
# htmlTable              2.1.0    2020-09-16 [2] CRAN (R 4.0.3)                           
# htmltools              0.5.0    2020-06-16 [1] CRAN (R 4.0.3)                           
# htmlwidgets            1.5.3    2020-12-10 [2] CRAN (R 4.0.3)                           
# httr                   1.4.2    2020-07-20 [1] CRAN (R 4.0.3)                           
# IRanges              * 2.24.0   2020-10-27 [1] Bioconductor                             
# iterators              1.0.13   2020-10-15 [2] CRAN (R 4.0.3)                           
# jaffelab             * 0.99.30  2020-11-02 [1] Github (LieberInstitute/jaffelab@42637ff)
# jpeg                   0.1-8.1  2019-10-24 [2] CRAN (R 4.0.3)                           
# jsonlite               1.7.1    2020-09-07 [1] CRAN (R 4.0.3)                           
# knitr                  1.30     2020-09-22 [1] CRAN (R 4.0.3)                           
# lattice                0.20-41  2020-04-02 [3] CRAN (R 4.0.3)                           
# latticeExtra           0.6-29   2019-12-19 [2] CRAN (R 4.0.3)                           
# lifecycle              0.2.0    2020-03-06 [1] CRAN (R 4.0.3)                           
# limma                  3.46.0   2020-10-27 [2] Bioconductor                             
# locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.0.3)                           
# lubridate              1.7.9.2  2020-11-13 [1] CRAN (R 4.0.3)                           
# magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.0.3)                           
# Matrix                 1.3-2    2021-01-06 [3] CRAN (R 4.0.3)                           
# MatrixGenerics       * 1.2.0    2020-10-27 [2] Bioconductor                             
# matrixStats          * 0.57.0   2020-09-25 [2] CRAN (R 4.0.3)                           
# memoise                1.1.0    2017-04-21 [2] CRAN (R 4.0.3)                           
# modelr                 0.1.8    2020-05-19 [2] CRAN (R 4.0.3)                           
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.0.3)                           
# nnet                   7.3-14   2020-04-26 [3] CRAN (R 4.0.3)                           
# openssl                1.4.3    2020-09-18 [1] CRAN (R 4.0.3)                           
# pheatmap             * 1.0.12   2019-01-04 [2] CRAN (R 4.0.3)                           
# pillar                 1.4.7    2020-11-20 [1] CRAN (R 4.0.3)                           
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.0.3)                           
# plyr                   1.8.6    2020-03-03 [2] CRAN (R 4.0.3)                           
# png                    0.1-7    2013-12-03 [2] CRAN (R 4.0.3)                           
# prettyunits            1.1.1    2020-01-24 [2] CRAN (R 4.0.3)                           
# progress               1.2.2    2019-05-16 [2] CRAN (R 4.0.3)                           
# ps                     1.4.0    2020-10-07 [1] CRAN (R 4.0.3)                           
# purrr                * 0.3.4    2020-04-17 [1] CRAN (R 4.0.3)                           
# qvalue                 2.22.0   2020-10-27 [2] Bioconductor                             
# R6                     2.5.0    2020-10-28 [1] CRAN (R 4.0.3)                           
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.0.3)                           
# rappdirs               0.3.1    2016-03-28 [2] CRAN (R 4.0.3)                           
# RColorBrewer         * 1.1-2    2014-12-07 [2] CRAN (R 4.0.3)                           
# Rcpp                   1.0.5    2020-07-06 [1] CRAN (R 4.0.3)                           
# RCurl                  1.98-1.2 2020-04-18 [2] CRAN (R 4.0.3)                           
# readr                * 1.4.0    2020-10-05 [2] CRAN (R 4.0.3)                           
# readxl                 1.3.1    2019-03-13 [2] CRAN (R 4.0.3)                           
# recount              * 1.16.1   2020-12-18 [2] Bioconductor                             
# rentrez                1.2.3    2020-11-10 [2] CRAN (R 4.0.3)                           
# reprex                 0.3.0    2019-05-16 [2] CRAN (R 4.0.3)                           
# reshape2             * 1.4.4    2020-04-09 [2] CRAN (R 4.0.3)                           
# rlang                * 0.4.9    2020-11-26 [1] CRAN (R 4.0.3)                           
# rngtools               1.5      2020-01-23 [2] CRAN (R 4.0.3)                           
# rpart                  4.1-15   2019-04-12 [3] CRAN (R 4.0.3)                           
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.0.3)                           
# Rsamtools              2.6.0    2020-10-27 [2] Bioconductor                             
# RSQLite                2.2.1    2020-09-30 [2] CRAN (R 4.0.3)                           
# rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.0.3)                           
# rtracklayer            1.50.0   2020-10-27 [2] Bioconductor                             
# rvest                  0.3.6    2020-07-25 [2] CRAN (R 4.0.3)                           
# S4Vectors            * 0.28.1   2020-12-09 [2] Bioconductor                             
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.0.3)                           
# scuttle              * 1.0.3    2020-11-23 [1] Bioconductor                             
# segmented              1.3-0    2020-10-27 [1] CRAN (R 4.0.3)                           
# sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.0.3)                           
# SingleCellExperiment * 1.12.0   2020-10-27 [2] Bioconductor                             
# sparseMatrixStats      1.2.0    2020-10-27 [2] Bioconductor                             
# stringi                1.5.3    2020-09-09 [1] CRAN (R 4.0.3)                           
# stringr              * 1.4.0    2019-02-10 [2] CRAN (R 4.0.3)                           
# SummarizedExperiment * 1.20.0   2020-10-27 [2] Bioconductor                             
# survival               3.2-7    2020-09-28 [3] CRAN (R 4.0.3)                           
# tibble               * 3.0.4    2020-10-12 [1] CRAN (R 4.0.3)                           
# tidyr                * 1.1.2    2020-08-27 [2] CRAN (R 4.0.3)                           
# tidyselect             1.1.0    2020-05-11 [2] CRAN (R 4.0.3)                           
# tidyverse            * 1.3.0    2019-11-21 [1] CRAN (R 4.0.3)                           
# utf8                   1.1.4    2018-05-24 [2] CRAN (R 4.0.3)                           
# VariantAnnotation      1.36.0   2020-10-27 [2] Bioconductor                             
# vctrs                  0.3.5    2020-11-17 [1] CRAN (R 4.0.3)                           
# withr                  2.3.0    2020-09-22 [1] CRAN (R 4.0.3)                           
# xfun                   0.19     2020-10-30 [1] CRAN (R 4.0.3)                           
# XML                    3.99-0.5 2020-07-23 [2] CRAN (R 4.0.3)                           
# xml2                   1.3.2    2020-04-23 [2] CRAN (R 4.0.3)                           
# XVector                0.30.0   2020-10-27 [2] Bioconductor                             
# zlibbioc               1.36.0   2020-10-27 [2] Bioconductor                             
# 
# [1] /users/lhuuki/R/4.0.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/library
