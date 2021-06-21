
library(SummarizedExperiment)
library(RColorBrewer)
library(tidyverse)
library(broom)
library(viridis)
library(DeconvoBuddies)
library(here)
library(sessioninfo)

## Load colors and plotting functions
source(here("main_colors.R"))
source(here("deconvolution","big_little_boxplot.R"))

cell_colors <- create_cell_colors(pallet = "classic")
region_colors <- list(Amygdala = "#FFFF1F",
                      sACC = "#8EB0F6")

## Load Bisque results & res data
load(here("deconvolution","data","est_prop_Bisque.Rdata"),verbose = TRUE)
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
load(here("deconvolution","data","sce_filtered.Rdata"), verbose = TRUE)

pd <- as.data.frame(colData(rse_gene))
pd2 <- pd[,c("RNum", "BrNum", "BrainRegion","Sex", "PrimaryDx", "Experiment")]
# Keep Bisque/All gene results
est_prop_bisque <- est_prop_bisque$all

long_prop <- est_prop_bisque$Est.prop.long %>%
  separate(sample, into = c("RNum", "Experiment"), extra = "merge") %>%
  left_join(pd2)

## Boxplots
boxplot_dx <- long_prop %>%
  ggplot(aes(x = cell_type, y = prop, fill = PrimaryDx)) +
  geom_boxplot() +
  labs(x = "Cell Type", y = "Proportion", fill ='Primary Dx') +
  theme_bw(base_size = 15)+
  scale_fill_manual(values = mdd_Dx_colors)+
  facet_wrap(~BrainRegion)

ggsave(boxplot_dx, filename = here("deconvolution", "plots", "bisque_cellType_boxplot_Dx.png"), width = 10)
ggsave(boxplot_dx, filename = here("deconvolution", "plots", "bisque_cellType_boxplot_Dx.pdf"), width = 10)

boxplot_region <- long_prop %>%
  ggplot(aes(x = cell_type, y = prop, fill = BrainRegion)) +
  geom_boxplot() +
  labs(x = "Cell Type", y = "Proportion", fill ='Brain Region') +
  scale_fill_manual(values = region_colors)+
  theme_bw(base_size = 15) 

ggsave(boxplot_region, filename = here("deconvolution", "plots", "bisque_cellType_boxplot_region.png"))
ggsave(boxplot_region, filename = here("deconvolution", "plots", "bisque_cellType_boxplot_region.pdf"))

boxplot_region <- long_prop %>%
  ggplot(aes(x = cell_type, y = prop, fill = BrainRegion)) +
  geom_boxplot() +
  labs(x = "Cell Type", y = "Proportion", fill ='Brain Region') +
  theme_bw(base_size = 15)+
  facet_wrap(~BrainRegion)

ggsave(boxplot_region, filename = here("deconvolution", "plots", "bisque_cellType_boxplot_region.png"))
ggsave(boxplot_region, filename = here("deconvolution", "plots", "bisque_cellType_boxplot_region.pdf"))


#### composition barplot ####
comp_barplot <- plot_composition_bar(long_prop, x_col = "BrainRegion") +
  scale_fill_manual(values = cell_colors)+
  theme_bw(base_size = 15)

ggsave(comp_barplot, filename = here("deconvolution","plots","bisque_composition_barplot.png"))
ggsave(comp_barplot, filename = here("deconvolution","plots","bisque_composition_barplot.pdf"))


# sgejobs::job_single('deconvo_plots', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript deconvo_plots.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value                                      
# version  R version 4.1.0 Patched (2021-05-18 r80330)
# os       CentOS Linux 7 (Core)                      
# system   x86_64, linux-gnu                          
# ui       X11                                        
# language (EN)                                       
# collate  en_US.UTF-8                                
# ctype    en_US.UTF-8                                
# tz       US/Eastern                                 
# date     2021-06-21                                 
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date       lib source                                  
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)                          
# backports              1.2.1    2020-12-09 [2] CRAN (R 4.1.0)                          
# beachmat               2.8.0    2021-05-19 [2] Bioconductor                            
# Biobase              * 2.52.0   2021-05-19 [2] Bioconductor                            
# BiocGenerics         * 0.38.0   2021-05-19 [2] Bioconductor                            
# BiocNeighbors          1.10.0   2021-05-19 [1] Bioconductor                            
# BiocParallel           1.26.0   2021-05-19 [2] Bioconductor                            
# BiocSingular           1.8.1    2021-06-08 [1] Bioconductor                            
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)                          
# bluster                1.2.1    2021-05-27 [1] Bioconductor                            
# broom                * 0.7.7    2021-06-13 [2] CRAN (R 4.1.0)                          
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.1.0)                          
# cli                    2.5.0    2021-04-26 [2] CRAN (R 4.1.0)                          
# cluster                2.1.2    2021-04-17 [3] CRAN (R 4.1.0)                          
# codetools              0.2-18   2020-11-04 [2] CRAN (R 4.1.0)                          
# colorout             * 1.2-2    2021-05-27 [1] Github (jalvesaq/colorout@79931fd)      
# colorspace             2.0-1    2021-05-04 [2] CRAN (R 4.1.0)                          
# crayon                 1.4.1    2021-02-08 [2] CRAN (R 4.1.0)                          
# DBI                    1.1.1    2021-01-15 [2] CRAN (R 4.1.0)                          
# dbplyr                 2.1.1    2021-04-06 [2] CRAN (R 4.1.0)                          
# DeconvoBuddies       * 0.99.0   2021-06-14 [1] Github (lahuuki/DeconvoBuddies@d889340) 
# DelayedArray           0.18.0   2021-05-19 [2] Bioconductor                            
# DelayedMatrixStats     1.14.0   2021-05-19 [2] Bioconductor                            
# digest                 0.6.27   2020-10-24 [2] CRAN (R 4.1.0)                          
# dplyr                * 1.0.7    2021-06-18 [2] CRAN (R 4.1.0)                          
# dqrng                  0.3.0    2021-05-01 [1] CRAN (R 4.1.0)                          
# edgeR                  3.34.0   2021-05-19 [2] Bioconductor                            
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)                          
# fansi                  0.5.0    2021-05-25 [2] CRAN (R 4.1.0)                          
# farver                 2.1.0    2021-02-28 [2] CRAN (R 4.1.0)                          
# forcats              * 0.5.1    2021-01-27 [2] CRAN (R 4.1.0)                          
# fs                     1.5.0    2020-07-31 [2] CRAN (R 4.1.0)                          
# generics               0.1.0    2020-10-31 [2] CRAN (R 4.1.0)                          
# GenomeInfoDb         * 1.28.0   2021-05-19 [2] Bioconductor                            
# GenomeInfoDbData       1.2.6    2021-05-11 [2] Bioconductor                            
# GenomicRanges        * 1.44.0   2021-05-19 [2] Bioconductor                            
# ggplot2              * 3.3.4    2021-06-16 [2] CRAN (R 4.1.0)                          
# glue                   1.4.2    2020-08-27 [2] CRAN (R 4.1.0)                          
# gridExtra              2.3      2017-09-09 [2] CRAN (R 4.1.0)                          
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)                          
# haven                  2.4.1    2021-04-23 [2] CRAN (R 4.1.0)                          
# here                 * 1.0.1    2020-12-13 [1] CRAN (R 4.1.0)                          
# hms                    1.1.0    2021-05-17 [2] CRAN (R 4.1.0)                          
# httr                   1.4.2    2020-07-20 [2] CRAN (R 4.1.0)                          
# igraph                 1.2.6    2020-10-06 [2] CRAN (R 4.1.0)                          
# IRanges              * 2.26.0   2021-05-19 [2] Bioconductor                            
# irlba                  2.3.3    2019-02-05 [2] CRAN (R 4.1.0)                          
# jsonlite               1.7.2    2020-12-09 [2] CRAN (R 4.1.0)                          
# labeling               0.4.2    2020-10-20 [2] CRAN (R 4.1.0)                          
# lattice                0.20-44  2021-05-02 [3] CRAN (R 4.1.0)                          
# lifecycle              1.0.0    2021-02-15 [2] CRAN (R 4.1.0)                          
# limma                  3.48.0   2021-05-19 [2] Bioconductor                            
# locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)                          
# lubridate              1.7.10   2021-02-26 [2] CRAN (R 4.1.0)                          
# magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.1.0)                          
# Matrix                 1.3-4    2021-06-01 [3] CRAN (R 4.1.0)                          
# MatrixGenerics       * 1.4.0    2021-05-19 [2] Bioconductor                            
# matrixStats          * 0.59.0   2021-06-01 [2] CRAN (R 4.1.0)                          
# metapod                1.0.0    2021-05-19 [1] Bioconductor                            
# modelr                 0.1.8    2020-05-19 [2] CRAN (R 4.1.0)                          
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)                          
# patchwork            * 1.1.1    2020-12-17 [1] CRAN (R 4.1.0)                          
# pillar                 1.6.1    2021-05-16 [2] CRAN (R 4.1.0)                          
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)                          
# pryr                   0.1.4    2018-02-18 [2] CRAN (R 4.1.0)                          
# ps                     1.6.0    2021-02-28 [2] CRAN (R 4.1.0)                          
# purrr                * 0.3.4    2020-04-17 [2] CRAN (R 4.1.0)                          
# R6                     2.5.0    2020-10-28 [2] CRAN (R 4.1.0)                          
# rafalib                1.0.0    2015-08-09 [1] CRAN (R 4.1.0)                          
# ragg                   1.1.3    2021-06-09 [1] CRAN (R 4.1.0)                          
# RColorBrewer         * 1.1-2    2014-12-07 [2] CRAN (R 4.1.0)                          
# Rcpp                   1.0.6    2021-01-15 [2] CRAN (R 4.1.0)                          
# RCurl                  1.98-1.3 2021-03-16 [2] CRAN (R 4.1.0)                          
# readr                * 1.4.0    2020-10-05 [2] CRAN (R 4.1.0)                          
# readxl                 1.3.1    2019-03-13 [2] CRAN (R 4.1.0)                          
# reprex                 2.0.0    2021-04-02 [2] CRAN (R 4.1.0)                          
# rlang                * 0.4.11   2021-04-30 [2] CRAN (R 4.1.0)                          
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)                          
# rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.1.0)                          
# rsvd                   1.0.5    2021-04-16 [1] CRAN (R 4.1.0)                          
# rvest                  1.0.0    2021-03-09 [2] CRAN (R 4.1.0)                          
# S4Vectors            * 0.30.0   2021-05-19 [2] Bioconductor                            
# ScaledMatrix           1.0.0    2021-05-19 [1] Bioconductor                            
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)                          
# scran                  1.20.1   2021-05-24 [1] Bioconductor                            
# scuttle                1.2.0    2021-05-19 [1] Bioconductor                            
# sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.1.0)                          
# sgejobs                0.99.1   2021-05-27 [1] Github (LieberInstitute/sgejobs@f5ab0ca)
# SingleCellExperiment   1.14.1   2021-05-21 [2] Bioconductor                            
# sparseMatrixStats      1.4.0    2021-05-19 [2] Bioconductor                            
# statmod                1.4.36   2021-05-10 [2] CRAN (R 4.1.0)                          
# stringi                1.6.2    2021-05-17 [2] CRAN (R 4.1.0)                          
# stringr              * 1.4.0    2019-02-10 [2] CRAN (R 4.1.0)                          
# SummarizedExperiment * 1.22.0   2021-05-19 [2] Bioconductor                            
# systemfonts            1.0.2    2021-05-11 [2] CRAN (R 4.1.0)                          
# textshaping            0.3.5    2021-06-09 [1] CRAN (R 4.1.0)                          
# tibble               * 3.1.2    2021-05-16 [2] CRAN (R 4.1.0)                          
# tidyr                * 1.1.3    2021-03-03 [2] CRAN (R 4.1.0)                          
# tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.1.0)                          
# tidyverse            * 1.3.1    2021-04-15 [2] CRAN (R 4.1.0)                          
# utf8                   1.2.1    2021-03-12 [2] CRAN (R 4.1.0)                          
# vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)                          
# viridis              * 0.6.1    2021-05-11 [2] CRAN (R 4.1.0)                          
# viridisLite          * 0.4.0    2021-04-13 [2] CRAN (R 4.1.0)                          
# withr                  2.4.2    2021-04-18 [2] CRAN (R 4.1.0)                          
# xml2                   1.3.2    2020-04-23 [2] CRAN (R 4.1.0)                          
# XVector                0.32.0   2021-05-19 [2] Bioconductor                            
# zlibbioc               1.38.0   2021-05-19 [2] Bioconductor                            
# 
# [1] /users/lhuuki/R/4.1
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/library
