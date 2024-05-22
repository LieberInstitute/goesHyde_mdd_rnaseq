
library("SummarizedExperiment")
library("SingleCellExperiment")
# library("jaffelab")
# library("xbioc")
library("BisqueRNA")
library("tidyverse")
library("DeconvoBuddies")
library("here")
library("sessioninfo")

#### Load Data ####
## Load rse_gene data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sce Data
load("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/deconvolution_bsp2/data/sce_pan.v2.Rdata", verbose = TRUE)

# 119 cell types
table(sce_pan$cellType)
table(sce_pan$cellType, sce_pan$cellType.Broad)
table(sce_pan$cellType.Broad)
table(is.na(sce_pan$cellType.Broad))
# Astro  Endo Macro Micro Mural Oligo   OPC Tcell Excit Inhib 
# 5596    31    32  3958   100 28165  4449    66  7617 20513 

#sACC & Amyg both have 5 donors (>4 min for Bisque)
table(sce_pan$donor, sce_pan$region)
#         amy dlpfc  hpc  nac sacc
# br5161 3269  4215 4352 2039 3159
# br5182    0     0    0 4209    0
# br5207    0  5294    0 4402    0
# br5212 3246  1693 3899 1753 3863
# br5276 2257     0    0 2522  588
# br5287    0     0 1853  669    0
# br5400 1736     0    0 4017 3930
# br5701 3498     0    0  281 3783

## 24 cell type for sACC (18 neurons)
table(droplevels(sce_pan[,sce_pan$region == "sacc"]$cellType))
# sacc_Astro_A sacc_Astro_B sacc_Excit_A sacc_Excit_B sacc_Excit_C sacc_Excit_D sacc_Excit_E sacc_Excit_F sacc_Excit_G 
# 747          160          856          575         1735          311          428          228           30 
# sacc_Inhib_A sacc_Inhib_B sacc_Inhib_C sacc_Inhib_D sacc_Inhib_E sacc_Inhib_F sacc_Inhib_G sacc_Inhib_H sacc_Inhib_I 
# 842          912          465          384          330          521          206          208           39 
# sacc_Inhib_J sacc_Inhib_K   sacc_Micro sacc_Oligo_A sacc_Oligo_B     sacc_OPC 
# 42           25          784         4389          195          911 

## 19 cell types for Amyg (11 neurons)
table(droplevels(sce_pan[,sce_pan$region == "amy"]$cellType))
# amy_Astro_A amy_Astro_B    amy_Endo amy_Excit_A amy_Excit_B amy_Excit_C amy_Inhib_A amy_Inhib_B amy_Inhib_C amy_Inhib_D 
# 1555          83          31         344          44          55         728         541         525         555 
# amy_Inhib_E amy_Inhib_F amy_Inhib_G amy_Inhib_H   amy_Micro   amy_Mural   amy_Oligo     amy_OPC   amy_Tcell 
# 414         216          86          52        1168          39        6080        1459          31


ratios <- get_mean_ratio2(sce_pan)
####


marker_stats %>%
  filter(rank_ratio <= 25, gene %in% rownames(rse_gene)) %>%
  count(cellType.target)

# # A tibble: 10 × 2
# cellType.target     n
# <fct>           <int>
#   1 Astro              25
# 2 Endo               25
# 3 Macro              21
# 4 Micro              25
# 5 Mural              24
# 6 Oligo              25
# 7 OPC                23
# 8 Tcell              14
# 9 Excit              25
# 10 Inhib              24

marker_stats2 <- marker_stats %>%
  filter(gene %in% rownames(rse_gene)) %>%
  filter(cellType.target != "Macro" | (cellType.target == "Macro" & rank_ratio <= 23)) %>%
  group_by(cellType.target) %>%
  mutate(rank_ratio = row_number())

marker_stats2%>%
  filter(rank_ratio <= 25) %>% 
  dplyr::count(cellType.target)

# # Groups:   cellType.target [10]
# cellType.target     n
# <fct>           <int>
#   1 Astro              25
# 2 Endo               25
# 3 Macro              18
# ...

marker_genes <- marker_stats2 %>%
  filter(rank_ratio <= 25) %>%
  pull(gene)

length(marker_genes)
# [1] 243

## save table of marker genes
## marker stats
rd <- as.data.frame(rowData(rse_gene)[marker_genes,]) %>% select(ensemblID, gencodeID)

marker_table <- marker_stats2 %>% filter(rank_ratio <= 25) %>%
  select(ensemblID = gene, `Cell Type` = cellType.target, Symbol) %>%
  left_join(rd) %>%
  select(gencodeID, `Cell Type`, Symbol)

write_csv(marker_table, file = here("deconvolution","data","MDDseq-deconvolution-markers.csv"))

#### create expression set ####
exp_set_bulk <- ExpressionSet(assayData = assays(rse_gene)$counts,
                              phenoData=AnnotatedDataFrame(
                                as.data.frame(colData(rse_gene))[c("BrNum")]))

exp_set_sce <- ExpressionSet(assayData = as.matrix(assays(sce_pan)$counts),
                             phenoData=AnnotatedDataFrame(
                               as.data.frame(colData(sce_pan))[c("cellType.Broad", "cellType", "uniqueID","donor")]))


exp_set_sce <- exp_set_sce[marker_genes,]
zero_cell_filter <- colSums(exprs(exp_set_sce)) != 0
message("Exclude ",sum(!zero_cell_filter), " cells")
# Exclude 2 cells

exp_set_sce <- exp_set_sce[,zero_cell_filter]

  #### Run Bisque ####
est_prop <- ReferenceBasedDecomposition(bulk.eset = exp_set_bulk,
                                        sc.eset = exp_set_sce,
                                        cell.types = "cellType.Broad",
                                        subject.names = "donor",
                                        use.overlap = FALSE)

est_prop$bulk.props <- t(est_prop$bulk.props)

est_prop$Est.prop.long <- est_prop$bulk.props %>%
  as.data.frame() %>%
  rownames_to_column("Sample")%>%
  pivot_longer(!Sample, names_to = "cell_type", values_to = "prop")

# est_prop$ilr <- ilr(est_prop$bulk.props)
# colnames(est_prop$ilr) <- paste0("ilr_",1:ncol(est_prop$ilr))

## Add long data and save
round(colMeans(est_prop$bulk.props), 3)
# Astro  Endo Macro Micro Mural Oligo   OPC Tcell Excit Inhib 
# 0.066 0.003 0.005 0.046 0.010 0.351 0.062 0.010 0.095 0.352 

est_prop_bisque <- est_prop
save(est_prop_bisque, file = here("deconvolution","data","est_prop_Bisque.Rdata"))

# sgejobs::job_single('deconvo_Bisque', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript deconvo_Bisque.R")
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
# date     2021-08-02                                 
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version    date       lib source                                   
# AnnotationDbi        * 1.54.1     2021-06-08 [2] Bioconductor                             
# assertthat             0.2.1      2019-03-21 [2] CRAN (R 4.1.0)                           
# backports              1.2.1      2020-12-09 [2] CRAN (R 4.1.0)                           
# Biobase              * 2.52.0     2021-05-19 [2] Bioconductor                             
# BiocGenerics         * 0.38.0     2021-05-19 [2] Bioconductor                             
# BiocManager            1.30.16    2021-06-15 [2] CRAN (R 4.1.0)                           
# Biostrings             2.60.1     2021-06-06 [2] Bioconductor                             
# BisqueRNA            * 1.0.5      2021-05-23 [1] CRAN (R 4.1.0)                           
# bit                    4.0.4      2020-08-04 [2] CRAN (R 4.1.0)                           
# bit64                  4.0.5      2020-08-30 [2] CRAN (R 4.1.0)                           
# bitops                 1.0-7      2021-04-24 [2] CRAN (R 4.1.0)                           
# blob                   1.2.2      2021-07-23 [2] CRAN (R 4.1.0)                           
# broom                  0.7.9      2021-07-27 [2] CRAN (R 4.1.0)                           
# cachem                 1.0.5      2021-05-15 [2] CRAN (R 4.1.0)                           
# cellranger             1.1.0      2016-07-27 [2] CRAN (R 4.1.0)                           
# checkmate              2.0.0      2020-02-06 [2] CRAN (R 4.1.0)                           
# cli                    3.0.1      2021-07-17 [2] CRAN (R 4.1.0)                           
# codetools              0.2-18     2020-11-04 [2] CRAN (R 4.1.0)                           
# colorout             * 1.2-2      2021-05-27 [1] Github (jalvesaq/colorout@79931fd)       
# colorspace             2.0-2      2021-06-24 [2] CRAN (R 4.1.0)                           
# crayon                 1.4.1      2021-02-08 [2] CRAN (R 4.1.0)                           
# DBI                    1.1.1      2021-01-15 [2] CRAN (R 4.1.0)                           
# dbplyr                 2.1.1      2021-04-06 [2] CRAN (R 4.1.0)                           
# DelayedArray           0.18.0     2021-05-19 [2] Bioconductor                             
# digest                 0.6.27     2020-10-24 [2] CRAN (R 4.1.0)                           
# dplyr                * 1.0.7      2021-06-18 [2] CRAN (R 4.1.0)                           
# ellipsis               0.3.2      2021-04-29 [2] CRAN (R 4.1.0)                           
# fansi                  0.5.0      2021-05-25 [2] CRAN (R 4.1.0)                           
# fastmap                1.1.0      2021-01-25 [2] CRAN (R 4.1.0)                           
# forcats              * 0.5.1      2021-01-27 [2] CRAN (R 4.1.0)                           
# fs                     1.5.0      2020-07-31 [2] CRAN (R 4.1.0)                           
# gargle                 1.2.0      2021-07-02 [2] CRAN (R 4.1.0)                           
# generics               0.1.0      2020-10-31 [2] CRAN (R 4.1.0)                           
# GenomeInfoDb         * 1.28.1     2021-07-01 [2] Bioconductor                             
# GenomeInfoDbData       1.2.6      2021-05-11 [2] Bioconductor                             
# GenomicRanges        * 1.44.0     2021-05-19 [2] Bioconductor                             
# ggplot2              * 3.3.5      2021-06-25 [2] CRAN (R 4.1.0)                           
# glue                   1.4.2      2020-08-27 [2] CRAN (R 4.1.0)                           
# googledrive            2.0.0      2021-07-08 [2] CRAN (R 4.1.0)                           
# gtable                 0.3.0      2019-03-25 [2] CRAN (R 4.1.0)                           
# haven                  2.4.1      2021-04-23 [2] CRAN (R 4.1.0)                           
# here                 * 1.0.1      2020-12-13 [1] CRAN (R 4.1.0)                           
# hms                    1.1.0      2021-05-17 [2] CRAN (R 4.1.0)                           
# httr                   1.4.2      2020-07-20 [2] CRAN (R 4.1.0)                           
# IRanges              * 2.26.0     2021-05-19 [2] Bioconductor                             
# jaffelab             * 0.99.31    2021-05-27 [1] Github (LieberInstitute/jaffelab@2cbd55a)
# jsonlite               1.7.2      2020-12-09 [2] CRAN (R 4.1.0)                           
# KEGGREST               1.32.0     2021-05-19 [2] Bioconductor                             
# lattice                0.20-44    2021-05-02 [3] CRAN (R 4.1.0)                           
# lifecycle              1.0.0      2021-02-15 [2] CRAN (R 4.1.0)                           
# limma                  3.48.1     2021-06-24 [2] Bioconductor                             
# limSolve               1.5.6      2019-11-14 [1] CRAN (R 4.1.0)                           
# lpSolve                5.6.15     2020-01-24 [1] CRAN (R 4.1.0)                           
# lubridate              1.7.10     2021-02-26 [2] CRAN (R 4.1.0)                           
# magrittr               2.0.1      2020-11-17 [2] CRAN (R 4.1.0)                           
# MASS                   7.3-54     2021-05-03 [3] CRAN (R 4.1.0)                           
# Matrix                 1.3-4      2021-06-01 [3] CRAN (R 4.1.0)                           
# MatrixGenerics       * 1.4.0      2021-05-19 [2] Bioconductor                             
# matrixStats          * 0.60.0     2021-07-26 [2] CRAN (R 4.1.0)                           
# memoise                2.0.0      2021-01-26 [2] CRAN (R 4.1.0)                           
# modelr                 0.1.8      2020-05-19 [2] CRAN (R 4.1.0)                           
# munsell                0.5.0      2018-06-12 [2] CRAN (R 4.1.0)                           
# pillar                 1.6.2      2021-07-29 [2] CRAN (R 4.1.0)                           
# pkgconfig              2.0.3      2019-09-22 [2] CRAN (R 4.1.0)                           
# pkgmaker               0.32.2.900 2021-06-24 [1] Github (renozao/pkgmaker@4818513)        
# plyr                   1.8.6      2020-03-03 [2] CRAN (R 4.1.0)                           
# png                    0.1-7      2013-12-03 [2] CRAN (R 4.1.0)                           
# purrr                * 0.3.4      2020-04-17 [2] CRAN (R 4.1.0)                           
# quadprog               1.5-8      2019-11-20 [2] CRAN (R 4.1.0)                           
# R6                     2.5.0      2020-10-28 [2] CRAN (R 4.1.0)                           
# rafalib              * 1.0.0      2015-08-09 [1] CRAN (R 4.1.0)                           
# RColorBrewer           1.1-2      2014-12-07 [2] CRAN (R 4.1.0)                           
# Rcpp                   1.0.7      2021-07-07 [2] CRAN (R 4.1.0)                           
# RCurl                  1.98-1.3   2021-03-16 [2] CRAN (R 4.1.0)                           
# readr                * 2.0.0      2021-07-20 [2] CRAN (R 4.1.0)                           
# readxl                 1.3.1      2019-03-13 [2] CRAN (R 4.1.0)                           
# registry               0.5-1      2019-03-05 [2] CRAN (R 4.1.0)                           
# reprex                 2.0.0      2021-04-02 [2] CRAN (R 4.1.0)                           
# reshape2               1.4.4      2020-04-09 [2] CRAN (R 4.1.0)                           
# rlang                  0.4.11     2021-04-30 [2] CRAN (R 4.1.0)                           
# rprojroot              2.0.2      2020-11-15 [2] CRAN (R 4.1.0)                           
# RSQLite                2.2.7      2021-04-22 [2] CRAN (R 4.1.0)                           
# rstudioapi             0.13       2020-11-12 [2] CRAN (R 4.1.0)                           
# rvest                  1.0.1      2021-07-26 [2] CRAN (R 4.1.0)                           
# S4Vectors            * 0.30.0     2021-05-19 [2] Bioconductor                             
# scales                 1.1.1      2020-05-11 [2] CRAN (R 4.1.0)                           
# segmented              1.3-4      2021-04-22 [1] CRAN (R 4.1.0)                           
# sessioninfo          * 1.1.1      2018-11-05 [2] CRAN (R 4.1.0)                           
# SingleCellExperiment * 1.14.1     2021-05-21 [2] Bioconductor                             
# stringi                1.7.3      2021-07-16 [2] CRAN (R 4.1.0)                           
# stringr              * 1.4.0      2019-02-10 [2] CRAN (R 4.1.0)                           
# SummarizedExperiment * 1.22.0     2021-05-19 [2] Bioconductor                             
# tibble               * 3.1.3      2021-07-23 [2] CRAN (R 4.1.0)                           
# tidyr                * 1.1.3      2021-03-03 [2] CRAN (R 4.1.0)                           
# tidyselect             1.1.1      2021-04-30 [2] CRAN (R 4.1.0)                           
# tidyverse            * 1.3.1      2021-04-15 [2] CRAN (R 4.1.0)                           
# tzdb                   0.1.2      2021-07-20 [2] CRAN (R 4.1.0)                           
# utf8                   1.2.2      2021-07-24 [2] CRAN (R 4.1.0)                           
# vctrs                  0.3.8      2021-04-29 [2] CRAN (R 4.1.0)                           
# vroom                  1.5.3      2021-07-14 [2] CRAN (R 4.1.0)                           
# withr                  2.4.2      2021-04-18 [2] CRAN (R 4.1.0)                           
# xbioc                * 0.1.19     2021-06-24 [1] Github (renozao/xbioc@1354168)           
# xml2                   1.3.2      2020-04-23 [2] CRAN (R 4.1.0)                           
# xtable                 1.8-4      2019-04-21 [2] CRAN (R 4.1.0)                           
# XVector                0.32.0     2021-05-19 [2] Bioconductor                             
# zlibbioc               1.38.0     2021-05-19 [2] Bioconductor                             
# 
# [1] /users/lhuuki/R/4.1
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/library


