library(SingleCellExperiment)
library(jaffelab)
library(here)
library(GenomicFeatures)
library(sessioninfo)

#### Load and filter data ####
## Load rse_gene data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sACC data
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose = TRUE)
sce.sacc$uniqueID <- paste0(sce.sacc$donor, "_", sce.sacc$Barcode)
colnames(sce.sacc) <- sce.sacc$uniqueID

sce.sacc <- sce.sacc[,sce.sacc$cellType != "Ambig.lowNtrxts",]
sce.sacc$cellType <- droplevels(sce.sacc$cellType)
## Add cellType.broad
sce.sacc$cellType.Broad <- ss(as.character(sce.sacc$cellType), "\\.", 1)
sce.sacc$cellType.Broad <- factor(sce.sacc$cellType.Broad)
## Match rownames
rownames(sce.sacc) <- rowData(sce.sacc)$ID
table(rownames(sce.sacc) %in% rownames(rse_gene))
# FALSE  TRUE 
# 15021 18517 

## filter for common between sc and bulk data
common_genes <- rownames(rse_gene)[rownames(rse_gene) %in% rownames(sce.sacc)]
sce.sacc <- sce.sacc[common_genes, ]

# Remove 0 genes across all nuclei
sce.sacc <- sce.sacc[!rowSums(assay(sce.sacc, "counts"))==0, ] 
dim(sce.sacc)
# [1] 17785  7004

## Amyg Data
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose = TRUE)
sce.amy$uniqueID <- paste0(sce.amy$donor, "_", sce.amy$Barcode)
colnames(sce.amy) <- sce.amy$uniqueID
## make cellType cols the same
sce.amy$cellType = sce.amy$cellType.split
sce.amy <- sce.amy[,sce.amy$cellType != "Ambig.lowNtrxts",]
sce.amy$cellType <- droplevels(sce.amy$cellType)
## Add cellType.broad
sce.amy$cellType.Broad <- ss(as.character(sce.amy$cellType), "\\.", 1)
sce.amy$cellType.Broad <- factor(sce.amy$cellType.Broad)
## Match rownames
rownames(sce.amy) <- rowData(sce.amy)$ID
table(rownames(sce.amy) %in% rownames(rse_gene))
# FALSE  TRUE 
# 15021 18517 

## filter for common between sc and bulk data
common_genes <- rownames(rse_gene)[rownames(rse_gene) %in% rownames(sce.amy)]
sce.amy <- sce.amy[common_genes, ]

# Remove 0 genes across all nuclei
sce.amy <- sce.amy[!rowSums(assay(sce.amy, "counts"))==0, ] 
dim(sce.amy)
# [1] 17792  6582

#### Add bp length for RPKM later #### 
rd.sacc <- rowData(sce.sacc)
rd.amy <- rowData(sce.amy)
## Import gnomic features
txdb <- makeTxDbFromGFF("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/GRCh38-3.0.0_premrna/genes/genes.gtf")
g <- genes(txdb)
e <- exonsBy(txdb, by = "gene")

if (!all(names(e) %in% names(g))) {
  warning("Dropping exons with gene ids not present in the gene list")
  e <- e[names(e) %in% names(g)]
}
e2 <- disjoin(e)
g$bp_length <- sum(width(e2))
summary(g$bp_length)

## Add data to sacc
g_sacc <- g[rownames(sce.sacc),]
summary(g_sacc$bp_length)
table(g_sacc$bp_length > 205012)

rowRanges(sce.sacc) <- g_sacc
all(rownames(rd.sacc) == rownames(sce.sacc))
rowData(sce.sacc)$Symbol <- rd.sacc$Symbol

## Add data to amy data
g_amy <- g[rownames(sce.amy),]
summary(g_amy$bp_length)
table(g_amy$bp_length > 205012)

rowRanges(sce.amy) <- g_amy
all(rownames(rd.amy) == rownames(sce.amy))
rowData(sce.amy)$Symbol <- rd.amy$Symbol

## Save filtered sce object
save(sce.sacc, file = here("deconvolution","data","sce.sacc_filtered.Rdata"))
## rename amy to amyg to be consistent with rest of mdd scripts
sce.amyg <- sce.amy
save(sce.amyg, file = here("deconvolution","data","sce.amyg_filtered.Rdata"))


# sgejobs::job_single('sce_data_prep', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript 1_sce_data_prep.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# [1] "2020-12-21 11:41:42 EST"
# user  system elapsed 
# 203.717   6.376 215.481 
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
# date     2020-12-21                                 
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date       lib source                                   
# AnnotationDbi        * 1.52.0   2020-10-27 [2] Bioconductor                             
# askpass                1.1      2019-01-13 [2] CRAN (R 4.0.3)                           
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.0.3)                           
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
# cli                    2.2.0    2020-11-20 [1] CRAN (R 4.0.3)                           
# crayon                 1.3.4    2017-09-16 [2] CRAN (R 4.0.3)                           
# curl                   4.3      2019-12-02 [2] CRAN (R 4.0.3)                           
# DBI                    1.1.0    2019-12-15 [2] CRAN (R 4.0.3)                           
# dbplyr                 2.0.0    2020-11-03 [2] CRAN (R 4.0.3)                           
# DelayedArray           0.16.0   2020-10-27 [2] Bioconductor                             
# digest                 0.6.27   2020-10-24 [1] CRAN (R 4.0.3)                           
# dplyr                  1.0.2    2020-08-18 [1] CRAN (R 4.0.3)                           
# ellipsis               0.3.1    2020-05-15 [2] CRAN (R 4.0.3)                           
# fansi                  0.4.1    2020-01-08 [2] CRAN (R 4.0.3)                           
# generics               0.1.0    2020-10-31 [2] CRAN (R 4.0.3)                           
# GenomeInfoDb         * 1.26.2   2020-12-08 [2] Bioconductor                             
# GenomeInfoDbData       1.2.4    2020-11-30 [2] Bioconductor                             
# GenomicAlignments      1.26.0   2020-10-27 [2] Bioconductor                             
# GenomicFeatures      * 1.42.1   2020-11-12 [2] Bioconductor                             
# GenomicRanges        * 1.42.0   2020-10-27 [2] Bioconductor                             
# glue                   1.4.2    2020-08-27 [1] CRAN (R 4.0.3)                           
# googledrive            1.0.1    2020-05-05 [1] CRAN (R 4.0.3)                           
# here                 * 1.0.0    2020-11-15 [1] CRAN (R 4.0.3)                           
# hms                    0.5.3    2020-01-08 [2] CRAN (R 4.0.3)                           
# httr                   1.4.2    2020-07-20 [1] CRAN (R 4.0.3)                           
# IRanges              * 2.24.0   2020-10-27 [1] Bioconductor                             
# jaffelab             * 0.99.30  2020-11-02 [1] Github (LieberInstitute/jaffelab@42637ff)
# lattice                0.20-41  2020-04-02 [3] CRAN (R 4.0.3)                           
# lifecycle              0.2.0    2020-03-06 [1] CRAN (R 4.0.3)                           
# limma                  3.46.0   2020-10-27 [2] Bioconductor                             
# magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.0.3)                           
# Matrix                 1.2-18   2019-11-27 [3] CRAN (R 4.0.3)                           
# MatrixGenerics       * 1.2.0    2020-10-27 [2] Bioconductor                             
# matrixStats          * 0.57.0   2020-09-25 [2] CRAN (R 4.0.3)                           
# memoise                1.1.0    2017-04-21 [2] CRAN (R 4.0.3)                           
# openssl                1.4.3    2020-09-18 [1] CRAN (R 4.0.3)                           
# pillar                 1.4.7    2020-11-20 [1] CRAN (R 4.0.3)                           
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.0.3)                           
# prettyunits            1.1.1    2020-01-24 [2] CRAN (R 4.0.3)                           
# progress               1.2.2    2019-05-16 [2] CRAN (R 4.0.3)                           
# purrr                  0.3.4    2020-04-17 [1] CRAN (R 4.0.3)                           
# R6                     2.5.0    2020-10-28 [1] CRAN (R 4.0.3)                           
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.0.3)                           
# rappdirs               0.3.1    2016-03-28 [2] CRAN (R 4.0.3)                           
# RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.0.3)                           
# Rcpp                   1.0.5    2020-07-06 [1] CRAN (R 4.0.3)                           
# RCurl                  1.98-1.2 2020-04-18 [2] CRAN (R 4.0.3)                           
# rlang                  0.4.9    2020-11-26 [1] CRAN (R 4.0.3)                           
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.0.3)                           
# Rsamtools              2.6.0    2020-10-27 [2] Bioconductor                             
# RSQLite                2.2.1    2020-09-30 [2] CRAN (R 4.0.3)                           
# rtracklayer            1.50.0   2020-10-27 [2] Bioconductor                             
# S4Vectors            * 0.28.1   2020-12-09 [2] Bioconductor                             
# segmented              1.3-0    2020-10-27 [1] CRAN (R 4.0.3)                           
# sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.0.3)                           
# SingleCellExperiment * 1.12.0   2020-10-27 [2] Bioconductor                             
# stringi                1.5.3    2020-09-09 [1] CRAN (R 4.0.3)                           
# stringr                1.4.0    2019-02-10 [2] CRAN (R 4.0.3)                           
# SummarizedExperiment * 1.20.0   2020-10-27 [2] Bioconductor                             
# tibble                 3.0.4    2020-10-12 [1] CRAN (R 4.0.3)                           
# tidyselect             1.1.0    2020-05-11 [2] CRAN (R 4.0.3)                           
# vctrs                  0.3.5    2020-11-17 [1] CRAN (R 4.0.3)                           
# withr                  2.3.0    2020-09-22 [1] CRAN (R 4.0.3)                           
# XML                    3.99-0.5 2020-07-23 [2] CRAN (R 4.0.3)                           
# xml2                   1.3.2    2020-04-23 [2] CRAN (R 4.0.3)                           
# XVector                0.30.0   2020-10-27 [2] Bioconductor                             
# zlibbioc               1.36.0   2020-10-27 [2] Bioconductor                             
# 
# [1] /users/lhuuki/R/4.0.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/library


