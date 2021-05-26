# qrsh -l bluejay,mem_free=50G,h_vmem=50G
library(SummarizedExperiment)
library(data.table)
library(edgeR)
library(recount)

# setwd("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline")

load("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/exprs_cutoff/rse_gene.Rdata")

gene_es <- as(rse_gene, "ExpressionSet")

assays(rse_gene)$counts <- getTPM(rse_gene, length_var = "Length", mapped_var = NULL)

tpm <- as(rse_gene, "ExpressionSet")

# lifted from https://rdrr.io/github/ctlab/phantasus/src/R/utils.R
write.gct <- function(es, file, gzip=FALSE) {
  if (gzip) {
    con <- gzfile(file)
  } else {
    con <- file(file)
  }
  open(con, open="w")
  writeLines("#1.3", con)
  ann.col <- ncol(pData(es))
  ann.row <- ncol(fData(es))
  writeLines(sprintf("%s\t%s\t%s\t%s", nrow(es), ncol(es), ann.row, ann.col), con)
  writeLines(paste0(c("ID", colnames(fData(es)), colnames(es)), collapse="\t"), con)
  
  ann.col.table <- t(as.matrix(pData(es)))
  ann.col.table <- cbind(matrix(rep(NA, ann.row*ann.col), nrow=ann.col), ann.col.table)
  write.table(ann.col.table, file=con, quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
  write.table(cbind(fData(es), exprs(es)), file=con, quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
  close(con)
}

write.gct(gene_es, here::here("predixcan_pipeline", "processed-data", "gene_counts.gct"))
write.gct(tpm, here::here("predixcan_pipeline", "processed-data", "tpm.gct"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

# > ## Reproducibility information
#   > print('Reproducibility information:')
# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2021-05-26 10:35:20 EDT"
# > proc.time()
# user  system elapsed 
# 160.436   8.541 923.038 
# > options(width = 120)
# > sessioninfo::session_info()
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
# date     2021-05-26                                 
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date       lib source                           
# abind                  1.4-5     2016-07-21 [2] CRAN (R 4.1.0)                   
# AnnotationDbi          1.54.0    2021-05-19 [2] Bioconductor                     
# assertthat             0.2.1     2019-03-21 [2] CRAN (R 4.1.0)                   
# backports              1.2.1     2020-12-09 [2] CRAN (R 4.1.0)                   
# base64enc              0.1-3     2015-07-28 [2] CRAN (R 4.1.0)                   
# Biobase              * 2.52.0    2021-05-19 [2] Bioconductor                     
# BiocFileCache          2.0.0     2021-05-19 [2] Bioconductor                     
# BiocGenerics         * 0.38.0    2021-05-19 [2] Bioconductor                     
# BiocIO                 1.2.0     2021-05-19 [2] Bioconductor                     
# BiocParallel           1.26.0    2021-05-19 [2] Bioconductor                     
# biomaRt                2.48.0    2021-05-19 [2] Bioconductor                     
# Biostrings             2.60.0    2021-05-19 [2] Bioconductor                     
# bit                    4.0.4     2020-08-04 [2] CRAN (R 4.1.0)                   
# bit64                  4.0.5     2020-08-30 [2] CRAN (R 4.1.0)                   
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.1.0)                   
# blob                   1.2.1     2020-01-20 [2] CRAN (R 4.1.0)                   
# boot                   1.3-28    2021-05-03 [3] CRAN (R 4.1.0)                   
# BSgenome               1.60.0    2021-05-19 [2] Bioconductor                     
# bumphunter             1.34.0    2021-05-19 [2] Bioconductor                     
# cachem                 1.0.5     2021-05-15 [2] CRAN (R 4.1.0)                   
# car                    3.0-10    2020-09-29 [2] CRAN (R 4.1.0)                   
# carData                3.0-4     2020-05-22 [2] CRAN (R 4.1.0)                   
# cellranger             1.1.0     2016-07-27 [2] CRAN (R 4.1.0)                   
# checkmate              2.0.0     2020-02-06 [2] CRAN (R 4.1.0)                   
# class                  7.3-19    2021-05-03 [3] CRAN (R 4.1.0)                   
# cli                    2.5.0     2021-04-26 [2] CRAN (R 4.1.0)                   
# cluster                2.1.2     2021-04-17 [3] CRAN (R 4.1.0)                   
# codetools              0.2-18    2020-11-04 [2] CRAN (R 4.1.0)                   
# colorspace             2.0-1     2021-05-04 [2] CRAN (R 4.1.0)                   
# crayon                 1.4.1     2021-02-08 [2] CRAN (R 4.1.0)                   
# curl                   4.3.1     2021-04-30 [2] CRAN (R 4.1.0)                   
# data.table           * 1.14.0    2021-02-21 [2] CRAN (R 4.1.0)                   
# DBI                    1.1.1     2021-01-15 [2] CRAN (R 4.1.0)                   
# dbplyr                 2.1.1     2021-04-06 [2] CRAN (R 4.1.0)                   
# DelayedArray           0.18.0    2021-05-19 [2] Bioconductor                     
# DEoptimR               1.0-9     2021-05-24 [2] CRAN (R 4.1.0)                   
# derfinder              1.26.0    2021-05-19 [2] Bioconductor                     
# derfinderHelper        1.26.0    2021-05-19 [2] Bioconductor                     
# destiny              * 3.1.1     2021-05-25 [1] Github (theislab/destiny@28307e9)
# digest                 0.6.27    2020-10-24 [2] CRAN (R 4.1.0)                   
# doRNG                  1.8.2     2020-01-27 [2] CRAN (R 4.1.0)                   
# downloader             0.4       2015-07-09 [2] CRAN (R 4.1.0)                   
# dplyr                  1.0.6     2021-05-05 [2] CRAN (R 4.1.0)                   
# e1071                  1.7-7     2021-05-23 [2] CRAN (R 4.1.0)                   
# edgeR                * 3.34.0    2021-05-19 [2] Bioconductor                     
# ellipsis               0.3.2     2021-04-29 [2] CRAN (R 4.1.0)                   
# fansi                  0.5.0     2021-05-25 [2] CRAN (R 4.1.0)                   
# fastmap                1.1.0     2021-01-25 [2] CRAN (R 4.1.0)                   
# filelock               1.0.2     2018-10-05 [2] CRAN (R 4.1.0)                   
# forcats                0.5.1     2021-01-27 [2] CRAN (R 4.1.0)                   
# foreach                1.5.1     2020-10-15 [2] CRAN (R 4.1.0)                   
# foreign                0.8-81    2020-12-22 [3] CRAN (R 4.1.0)                   
# Formula                1.2-4     2020-10-16 [2] CRAN (R 4.1.0)                   
# generics               0.1.0     2020-10-31 [2] CRAN (R 4.1.0)                   
# GenomeInfoDb         * 1.28.0    2021-05-19 [2] Bioconductor                     
# GenomeInfoDbData       1.2.6     2021-05-11 [2] Bioconductor                     
# GenomicAlignments      1.28.0    2021-05-19 [2] Bioconductor                     
# GenomicFeatures        1.44.0    2021-05-19 [2] Bioconductor                     
# GenomicFiles           1.28.0    2021-05-19 [2] Bioconductor                     
# GenomicRanges        * 1.44.0    2021-05-19 [2] Bioconductor                     
# GEOquery               2.60.0    2021-05-19 [2] Bioconductor                     
# ggplot.multistats      1.0.0     2019-10-28 [1] CRAN (R 4.1.0)                   
# ggplot2                3.3.3     2020-12-30 [2] CRAN (R 4.1.0)                   
# ggthemes               4.2.4     2021-01-20 [1] CRAN (R 4.1.0)                   
# glue                   1.4.2     2020-08-27 [2] CRAN (R 4.1.0)                   
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.1.0)                   
# gtable                 0.3.0     2019-03-25 [2] CRAN (R 4.1.0)                   
# haven                  2.4.1     2021-04-23 [2] CRAN (R 4.1.0)                   
# here                   1.0.1     2020-12-13 [1] CRAN (R 4.1.0)                   
# hexbin                 1.28.2    2021-01-08 [2] CRAN (R 4.1.0)                   
# Hmisc                  4.5-0     2021-02-28 [2] CRAN (R 4.1.0)                   
# hms                    1.1.0     2021-05-17 [2] CRAN (R 4.1.0)                   
# htmlTable              2.2.1     2021-05-18 [2] CRAN (R 4.1.0)                   
# htmltools              0.5.1.1   2021-01-22 [2] CRAN (R 4.1.0)                   
# htmlwidgets            1.5.3     2020-12-10 [2] CRAN (R 4.1.0)                   
# httr                   1.4.2     2020-07-20 [2] CRAN (R 4.1.0)                   
# IRanges              * 2.26.0    2021-05-19 [2] Bioconductor                     
# irlba                  2.3.3     2019-02-05 [2] CRAN (R 4.1.0)                   
# iterators              1.0.13    2020-10-15 [2] CRAN (R 4.1.0)                   
# jpeg                   0.1-8.1   2019-10-24 [2] CRAN (R 4.1.0)                   
# jsonlite               1.7.2     2020-12-09 [2] CRAN (R 4.1.0)                   
# KEGGREST               1.32.0    2021-05-19 [2] Bioconductor                     
# knitr                  1.33      2021-04-24 [2] CRAN (R 4.1.0)                   
# laeken                 0.5.1     2020-02-05 [1] CRAN (R 4.1.0)                   
# lattice                0.20-44   2021-05-02 [3] CRAN (R 4.1.0)                   
# latticeExtra           0.6-29    2019-12-19 [2] CRAN (R 4.1.0)                   
# lifecycle              1.0.0     2021-02-15 [2] CRAN (R 4.1.0)                   
# limma                * 3.48.0    2021-05-19 [2] Bioconductor                     
# lmtest                 0.9-38    2020-09-09 [2] CRAN (R 4.1.0)                   
# locfit                 1.5-9.4   2020-03-25 [2] CRAN (R 4.1.0)                   
# magrittr               2.0.1     2020-11-17 [2] CRAN (R 4.1.0)                   
# MASS                   7.3-54    2021-05-03 [3] CRAN (R 4.1.0)                   
# Matrix                 1.3-3     2021-05-04 [3] CRAN (R 4.1.0)                   
# MatrixGenerics       * 1.4.0     2021-05-19 [2] Bioconductor                     
# matrixStats          * 0.58.0    2021-01-29 [2] CRAN (R 4.1.0)                   
# memoise                2.0.0     2021-01-26 [2] CRAN (R 4.1.0)                   
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.1.0)                   
# nnet                   7.3-16    2021-05-03 [3] CRAN (R 4.1.0)                   
# openxlsx               4.2.3     2020-10-27 [2] CRAN (R 4.1.0)                   
# pcaMethods             1.84.0    2021-05-19 [2] Bioconductor                     
# pillar                 1.6.1     2021-05-16 [2] CRAN (R 4.1.0)                   
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.1.0)                   
# plyr                   1.8.6     2020-03-03 [2] CRAN (R 4.1.0)                   
# png                    0.1-7     2013-12-03 [2] CRAN (R 4.1.0)                   
# prettyunits            1.1.1     2020-01-24 [2] CRAN (R 4.1.0)                   
# progress               1.2.2     2019-05-16 [2] CRAN (R 4.1.0)                   
# proxy                  0.4-25    2021-03-05 [2] CRAN (R 4.1.0)                   
# purrr                  0.3.4     2020-04-17 [2] CRAN (R 4.1.0)                   
# qvalue                 2.24.0    2021-05-19 [2] Bioconductor                     
# R6                     2.5.0     2020-10-28 [2] CRAN (R 4.1.0)                   
# ranger                 0.12.1    2020-01-10 [1] CRAN (R 4.1.0)                   
# rappdirs               0.3.3     2021-01-31 [2] CRAN (R 4.1.0)                   
# RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 4.1.0)                   
# Rcpp                   1.0.6     2021-01-15 [2] CRAN (R 4.1.0)                   
# RcppEigen              0.3.3.9.1 2020-12-17 [2] CRAN (R 4.1.0)                   
# RcppHNSW               0.3.0     2020-09-06 [1] CRAN (R 4.1.0)                   
# RCurl                  1.98-1.3  2021-03-16 [2] CRAN (R 4.1.0)                   
# readr                  1.4.0     2020-10-05 [2] CRAN (R 4.1.0)                   
# readxl                 1.3.1     2019-03-13 [2] CRAN (R 4.1.0)                   
# recount              * 1.18.0    2021-05-19 [2] Bioconductor                     
# rentrez                1.2.3     2020-11-10 [2] CRAN (R 4.1.0)                   
# reshape2               1.4.4     2020-04-09 [2] CRAN (R 4.1.0)                   
# restfulr               0.0.13    2017-08-06 [2] CRAN (R 4.1.0)                   
# rio                    0.5.26    2021-03-01 [2] CRAN (R 4.1.0)                   
# rjson                  0.2.20    2018-06-08 [2] CRAN (R 4.1.0)                   
# rlang                  0.4.11    2021-04-30 [2] CRAN (R 4.1.0)                   
# rngtools               1.5       2020-01-23 [2] CRAN (R 4.1.0)                   
# robustbase             0.93-7    2021-01-04 [2] CRAN (R 4.1.0)                   
# rpart                  4.1-15    2019-04-12 [3] CRAN (R 4.1.0)                   
# rprojroot              2.0.2     2020-11-15 [2] CRAN (R 4.1.0)                   
# Rsamtools              2.8.0     2021-05-19 [2] Bioconductor                     
# RSpectra               0.16-0    2019-12-01 [2] CRAN (R 4.1.0)                   
# RSQLite                2.2.7     2021-04-22 [2] CRAN (R 4.1.0)                   
# rstudioapi             0.13      2020-11-12 [2] CRAN (R 4.1.0)                   
# rtracklayer            1.52.0    2021-05-19 [2] Bioconductor                     
# S4Vectors            * 0.30.0    2021-05-19 [2] Bioconductor                     
# scales                 1.1.1     2020-05-11 [2] CRAN (R 4.1.0)                   
# scatterplot3d          0.3-41    2018-03-14 [1] CRAN (R 4.1.0)                   
# sessioninfo            1.1.1     2018-11-05 [2] CRAN (R 4.1.0)                   
# SingleCellExperiment   1.14.1    2021-05-21 [2] Bioconductor                     
# smoother               1.1       2015-04-16 [1] CRAN (R 4.1.0)                   
# sp                     1.4-5     2021-01-10 [2] CRAN (R 4.1.0)                   
# stringi                1.6.2     2021-05-17 [2] CRAN (R 4.1.0)                   
# stringr                1.4.0     2019-02-10 [2] CRAN (R 4.1.0)                   
# SummarizedExperiment * 1.22.0    2021-05-19 [2] Bioconductor                     
# survival               3.2-11    2021-04-26 [3] CRAN (R 4.1.0)                   
# tibble                 3.1.2     2021-05-16 [2] CRAN (R 4.1.0)                   
# tidyr                  1.1.3     2021-03-03 [2] CRAN (R 4.1.0)                   
# tidyselect             1.1.1     2021-04-30 [2] CRAN (R 4.1.0)                   
# TTR                    0.24.2    2020-09-01 [1] CRAN (R 4.1.0)                   
# utf8                   1.2.1     2021-03-12 [2] CRAN (R 4.1.0)                   
# VariantAnnotation      1.38.0    2021-05-19 [2] Bioconductor                     
# vcd                    1.4-8     2020-09-21 [1] CRAN (R 4.1.0)                   
# vctrs                  0.3.8     2021-04-29 [2] CRAN (R 4.1.0)                   
# VIM                    6.1.0     2021-01-19 [1] CRAN (R 4.1.0)                   
# withr                  2.4.2     2021-04-18 [2] CRAN (R 4.1.0)                   
# xfun                   0.23      2021-05-15 [2] CRAN (R 4.1.0)                   
# XML                    3.99-0.6  2021-03-16 [2] CRAN (R 4.1.0)                   
# xml2                   1.3.2     2020-04-23 [2] CRAN (R 4.1.0)                   
# xts                    0.12.1    2020-09-09 [2] CRAN (R 4.1.0)                   
# XVector                0.32.0    2021-05-19 [2] Bioconductor                     
# yaml                   2.2.1     2020-02-01 [2] CRAN (R 4.1.0)                   
# zip                    2.1.1     2020-08-27 [2] CRAN (R 4.1.0)                   
# zlibbioc               1.38.0    2021-05-19 [2] Bioconductor                     
# zoo                    1.8-9     2021-03-09 [2] CRAN (R 4.1.0)                   
# 
# [1] /users/aseyedia/R/4.1
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/library"