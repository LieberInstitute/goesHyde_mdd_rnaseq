# qrsh -l bluejay,mem_free=100G,h_vmem=100G
library(SummarizedExperiment)
library(data.table)
library(edgeR)
library(recount)

# setwd("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline")
# cd /dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline

load(here::here("exprs_cutoff", "rse_gene.Rdata"))

gene_es <- as(rse_gene, "ExpressionSet")

assays(rse_gene)$counts <- getTPM(rse_gene, length_var = "Length", mapped_var = NULL)

tpm <- as(rse_gene, "ExpressionSet")

vcf_samp <- readLines(here::here("predixcan_pipeline", "processed-data", "01_get_inv_quantile_norm", "vcf_samples.txt"))

# all TRUE
# colData(rse_gene)$genoSample[colData(rse_gene)$genoSample %in% vcf_samp]
# table(colData(rse_gene)$genoSample %in% vcf_samp)
# unique(data.table(colData(rse_gene)$genoSample, colData(rse_gene)$genoSample))
# data.table(unique(colData(rse_gene)$genoSample[colData(rse_gene)$genoSample %in% vcf_samp]), vcf_samp[vcf_samp %in% colData(rse_gene)$genoSample])

# new_df = df[vcf_samp, ]

# make one list of samples that have sample IDs in both the VCF and the expression and then reorder from there
# PLINK should be able to reorder okay
# you might have to make the vcf into a PLINK file before reordering but not necessary to copy the file

# might as well throw in the filtering at that point (maf) ask Leo about this

if(all(colData(rse_gene)$genoSample %in% vcf_samp)) {
  samp_lookup <- unique(data.table(colData(rse_gene)$genoSample, colData(rse_gene)$genoSample))
  colnames(samp_lookup) <- c("sample_id", "participant_id")
  write.table(samp_lookup, quote=FALSE, row.names = FALSE, file = here::here("predixcan_pipeline", "processed-data", "01_get_inv_quantile_norm", "samp_part_lookup.txt"))
}

# vcf_samp[vcf_samp %in% colData(rse_gene)$genoSample]

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

write.gct(gene_es, here::here("predixcan_pipeline", "processed-data", "01_get_inv_quantile_norm", "gene_counts.gct"))
write.gct(tpm, here::here("predixcan_pipeline", "processed-data", "01_get_inv_quantile_norm", "tpm.gct"))

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
# [1] "2021-06-02 13:15:14 EDT"
# > proc.time()
# user   system  elapsed 
# 30.480    2.618 1168.798 
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
# date     2021-06-02                                 
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date       lib source        
# AnnotationDbi          1.54.0   2021-05-19 [2] Bioconductor  
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# backports              1.2.1    2020-12-09 [2] CRAN (R 4.1.0)
# base64enc              0.1-3    2015-07-28 [2] CRAN (R 4.1.0)
# Biobase              * 2.52.0   2021-05-19 [2] Bioconductor  
# BiocFileCache          2.0.0    2021-05-19 [2] Bioconductor  
# BiocGenerics         * 0.38.0   2021-05-19 [2] Bioconductor  
# BiocIO                 1.2.0    2021-05-19 [2] Bioconductor  
# BiocParallel           1.26.0   2021-05-19 [2] Bioconductor  
# biomaRt                2.48.0   2021-05-19 [2] Bioconductor  
# Biostrings             2.60.0   2021-05-19 [2] Bioconductor  
# bit                    4.0.4    2020-08-04 [2] CRAN (R 4.1.0)
# bit64                  4.0.5    2020-08-30 [2] CRAN (R 4.1.0)
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# blob                   1.2.1    2020-01-20 [2] CRAN (R 4.1.0)
# BSgenome               1.60.0   2021-05-19 [2] Bioconductor  
# bumphunter             1.34.0   2021-05-19 [2] Bioconductor  
# cachem                 1.0.5    2021-05-15 [2] CRAN (R 4.1.0)
# checkmate              2.0.0    2020-02-06 [2] CRAN (R 4.1.0)
# cli                    2.5.0    2021-04-26 [2] CRAN (R 4.1.0)
# cluster                2.1.2    2021-04-17 [3] CRAN (R 4.1.0)
# codetools              0.2-18   2020-11-04 [2] CRAN (R 4.1.0)
# colorspace             2.0-1    2021-05-04 [2] CRAN (R 4.1.0)
# crayon                 1.4.1    2021-02-08 [2] CRAN (R 4.1.0)
# curl                   4.3.1    2021-04-30 [2] CRAN (R 4.1.0)
# data.table           * 1.14.0   2021-02-21 [2] CRAN (R 4.1.0)
# DBI                    1.1.1    2021-01-15 [2] CRAN (R 4.1.0)
# dbplyr                 2.1.1    2021-04-06 [2] CRAN (R 4.1.0)
# DelayedArray           0.18.0   2021-05-19 [2] Bioconductor  
# derfinder              1.26.0   2021-05-19 [2] Bioconductor  
# derfinderHelper        1.26.0   2021-05-19 [2] Bioconductor  
# digest                 0.6.27   2020-10-24 [2] CRAN (R 4.1.0)
# doRNG                  1.8.2    2020-01-27 [2] CRAN (R 4.1.0)
# downloader             0.4      2015-07-09 [2] CRAN (R 4.1.0)
# dplyr                  1.0.6    2021-05-05 [2] CRAN (R 4.1.0)
# edgeR                * 3.34.0   2021-05-19 [2] Bioconductor  
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  0.5.0    2021-05-25 [2] CRAN (R 4.1.0)
# fastmap                1.1.0    2021-01-25 [2] CRAN (R 4.1.0)
# filelock               1.0.2    2018-10-05 [2] CRAN (R 4.1.0)
# foreach                1.5.1    2020-10-15 [2] CRAN (R 4.1.0)
# foreign                0.8-81   2020-12-22 [3] CRAN (R 4.1.0)
# Formula                1.2-4    2020-10-16 [2] CRAN (R 4.1.0)
# generics               0.1.0    2020-10-31 [2] CRAN (R 4.1.0)
# GenomeInfoDb         * 1.28.0   2021-05-19 [2] Bioconductor  
# GenomeInfoDbData       1.2.6    2021-05-11 [2] Bioconductor  
# GenomicAlignments      1.28.0   2021-05-19 [2] Bioconductor  
# GenomicFeatures        1.44.0   2021-05-19 [2] Bioconductor  
# GenomicFiles           1.28.0   2021-05-19 [2] Bioconductor  
# GenomicRanges        * 1.44.0   2021-05-19 [2] Bioconductor  
# GEOquery               2.60.0   2021-05-19 [2] Bioconductor  
# ggplot2                3.3.3    2020-12-30 [2] CRAN (R 4.1.0)
# glue                   1.4.2    2020-08-27 [2] CRAN (R 4.1.0)
# gridExtra              2.3      2017-09-09 [2] CRAN (R 4.1.0)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# here                   1.0.1    2020-12-13 [1] CRAN (R 4.1.0)
# Hmisc                  4.5-0    2021-02-28 [2] CRAN (R 4.1.0)
# hms                    1.1.0    2021-05-17 [2] CRAN (R 4.1.0)
# htmlTable              2.2.1    2021-05-18 [2] CRAN (R 4.1.0)
# htmltools              0.5.1.1  2021-01-22 [2] CRAN (R 4.1.0)
# htmlwidgets            1.5.3    2020-12-10 [2] CRAN (R 4.1.0)
# httr                   1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
# IRanges              * 2.26.0   2021-05-19 [2] Bioconductor  
# iterators              1.0.13   2020-10-15 [2] CRAN (R 4.1.0)
# jpeg                   0.1-8.1  2019-10-24 [2] CRAN (R 4.1.0)
# jsonlite               1.7.2    2020-12-09 [2] CRAN (R 4.1.0)
# KEGGREST               1.32.0   2021-05-19 [2] Bioconductor  
# knitr                  1.33     2021-04-24 [2] CRAN (R 4.1.0)
# lattice                0.20-44  2021-05-02 [3] CRAN (R 4.1.0)
# latticeExtra           0.6-29   2019-12-19 [2] CRAN (R 4.1.0)
# lifecycle              1.0.0    2021-02-15 [2] CRAN (R 4.1.0)
# limma                * 3.48.0   2021-05-19 [2] Bioconductor  
# locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
# magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
# Matrix                 1.3-4    2021-06-01 [3] CRAN (R 4.1.0)
# MatrixGenerics       * 1.4.0    2021-05-19 [2] Bioconductor  
# matrixStats          * 0.59.0   2021-06-01 [2] CRAN (R 4.1.0)
# memoise                2.0.0    2021-01-26 [2] CRAN (R 4.1.0)
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# nnet                   7.3-16   2021-05-03 [3] CRAN (R 4.1.0)
# pillar                 1.6.1    2021-05-16 [2] CRAN (R 4.1.0)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# plyr                   1.8.6    2020-03-03 [2] CRAN (R 4.1.0)
# png                    0.1-7    2013-12-03 [2] CRAN (R 4.1.0)
# prettyunits            1.1.1    2020-01-24 [2] CRAN (R 4.1.0)
# progress               1.2.2    2019-05-16 [2] CRAN (R 4.1.0)
# purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# qvalue                 2.24.0   2021-05-19 [2] Bioconductor  
# R6                     2.5.0    2020-10-28 [2] CRAN (R 4.1.0)
# rappdirs               0.3.3    2021-01-31 [2] CRAN (R 4.1.0)
# RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
# Rcpp                   1.0.6    2021-01-15 [2] CRAN (R 4.1.0)
# RCurl                  1.98-1.3 2021-03-16 [2] CRAN (R 4.1.0)
# readr                  1.4.0    2020-10-05 [2] CRAN (R 4.1.0)
# recount              * 1.18.0   2021-05-19 [2] Bioconductor  
# rentrez                1.2.3    2020-11-10 [2] CRAN (R 4.1.0)
# reshape2               1.4.4    2020-04-09 [2] CRAN (R 4.1.0)
# restfulr               0.0.13   2017-08-06 [2] CRAN (R 4.1.0)
# rjson                  0.2.20   2018-06-08 [2] CRAN (R 4.1.0)
# rlang                  0.4.11   2021-04-30 [2] CRAN (R 4.1.0)
# rngtools               1.5      2020-01-23 [2] CRAN (R 4.1.0)
# rpart                  4.1-15   2019-04-12 [3] CRAN (R 4.1.0)
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
# Rsamtools              2.8.0    2021-05-19 [2] Bioconductor  
# RSQLite                2.2.7    2021-04-22 [2] CRAN (R 4.1.0)
# rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.1.0)
# rtracklayer            1.52.0   2021-05-19 [2] Bioconductor  
# S4Vectors            * 0.30.0   2021-05-19 [2] Bioconductor  
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
# sessioninfo            1.1.1    2018-11-05 [2] CRAN (R 4.1.0)
# stringi                1.6.2    2021-05-17 [2] CRAN (R 4.1.0)
# stringr                1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment * 1.22.0   2021-05-19 [2] Bioconductor  
# survival               3.2-11   2021-04-26 [3] CRAN (R 4.1.0)
# tibble                 3.1.2    2021-05-16 [2] CRAN (R 4.1.0)
# tidyr                  1.1.3    2021-03-03 [2] CRAN (R 4.1.0)
# tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
# utf8                   1.2.1    2021-03-12 [2] CRAN (R 4.1.0)
# VariantAnnotation      1.38.0   2021-05-19 [2] Bioconductor  
# vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
# withr                  2.4.2    2021-04-18 [2] CRAN (R 4.1.0)
# xfun                   0.23     2021-05-15 [2] CRAN (R 4.1.0)
# XML                    3.99-0.6 2021-03-16 [2] CRAN (R 4.1.0)
# xml2                   1.3.2    2020-04-23 [2] CRAN (R 4.1.0)
# XVector                0.32.0   2021-05-19 [2] Bioconductor  
# yaml                   2.2.1    2020-02-01 [2] CRAN (R 4.1.0)
# zlibbioc               1.38.0   2021-05-19 [2] Bioconductor  
# 
# [1] /users/aseyedia/R/4.1
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/library