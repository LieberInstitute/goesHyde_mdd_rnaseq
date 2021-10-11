#######
library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(WGCNA)
library(sva)
library(sessioninfo)
library(here)

#qrsh -l mem_free=100G,h_vmem=100G

## multithread
#setallowWGCNAThreads(8)

options(stringsAsFactors = FALSE)
#enableWGCNAThreads(8)


#dir.create("rdas", showWarnings = FALSE)

## load data
#load(here("exprs_cutoff","rse_gene.Rdata"))

load("../exprs_cutoff/rse_gene.Rdata", verbose = TRUE)
load("../data/degradation_rse_MDDseq_BiPSeq_BothRegions.Rdata", verbose = TRUE)

## checks
identical(colnames(rse_gene), colnames(cov_rse)) # TRUE
table(rse_gene$PrimaryDx)
# MDD Control Bipolar
# 463     387     245

table(cov_rse$PrimaryDx)

rse_gene <- rse_gene[, rse_gene$PrimaryDx %in% c("Control", "MDD")]
cov_rse <- cov_rse[, cov_rse$PrimaryDx %in% c("Control", "MDD")]
table(rse_gene$PrimaryDx)
# MDD Control Bipolar
# 463     387       0

rse_gene$Dx <- droplevels(rse_gene$PrimaryDx)
cov_rse$Dx <- droplevels(cov_rse$PrimaryDx)

table(rse_gene$Dx)
# MDD Control
# 463     387

## add ancestry
#no need -- already in new rse_object
#load("../genotype_data/goesHyde_bipolarMdd_Genotypes_n588_mds.rda", verbose = TRUE)
#class(mds)
#dim(mds)
#corner(mds)
#table(rse_gene$BrNum %in% rownames(mds))

# FALSE  TRUE
# 12   836
## 6 individual, 12 brain regions absent in mds file

#table(rse_gene$BrNum %in% rownames(mds), rse_gene$BrainRegion)
# Amygdala sACC
# FALSE        6    6
# TRUE       415  421

## keep samples with genotypes
#rse_gene <- rse_gene[, rse_gene$BrNum %in% rownames(mds)]
#cov_rse <- cov_rse[, cov_rse$BrNum %in% rownames(mds)]

dim(rse_gene)
# [1] 25212   846
addmargins(table(rse_gene$Dx, rse_gene$BrainRegion))
#         Amygdala sACC Sum
#  MDD          231  228 459
#  Control      187  200 387
#  Sum          418  428 846

#### reorders according to rse_gene$BrNum
#mds = mds[rse_gene$BrNum,1:5]
#dim(mds)
# [1] 836   5

#colnames(mds) = paste0("snpPC", 1:5)
#colData(rse_gene) = cbind(colData(rse_gene), mds)
#colData(cov_rse) = cbind(colData(cov_rse), mds)

###########
#compute RPKM
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')


##########
## model #
##########

### came from qSV_model_DE_analysis.R
# modJoint = model.matrix(~PrimaryDx + AgeDeath + Sex + mitoRate + rRNA_rate +
                            # totalAssignedGene + RIN, data = colData(rse_gene))


### removed ERCCsumLogErr term
### REPLACED ERCCsumLogErr term

modJoint = model.matrix(~Dx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 +
	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr,
	data=colData(rse_gene))


### counts from degrafation into log2 scale
degExprs = log2(assays(cov_rse)$count+1)
k = num.sv(degExprs, modJoint)
print(k) # 22
qSV_mat = prcomp(t(degExprs))$x[,1:k]

## join and move around region, dx and interaction for cleaning

colnames(modJoint)
# [1] "(Intercept)"               "DxControl"
# [3] "BrainRegionsACC"           "AgeDeath"
# [5] "SexM"                      "snpPC1"
# [7] "snpPC2"                    "snpPC3"
# [9] "mitoRate"                  "rRNA_rate"
#[11] "totalAssignedGene"         "RIN"
#[13] "ERCCsumLogErr"             "DxControl:BrainRegionsACC"

modQsva = cbind(modJoint[,c(1:4,14,5:13)], qSV_mat)

colnames(modQsva)

 #[1] "(Intercept)"               "DxControl"
 #[3] "BrainRegionsACC"           "AgeDeath"
 #[5] "DxControl:BrainRegionsACC" "SexM"
 #[7] "snpPC1"                    "snpPC2"
 #[9] "snpPC3"                    "mitoRate"
#[11] "rRNA_rate"                 "totalAssignedGene"
#[13] "RIN"                       "ERCCsumLogErr"
#[15] "PC1"                       "PC2"
#[17] "PC3"                       "PC4"
#[19] "PC5"                       "PC6"
#[21] "PC7"                       "PC8"
#[23] "PC9"                       "PC10"
#[25] "PC11"                      "PC12"
#[27] "PC13"                      "PC14"
#[29] "PC15"                      "PC16"
#[31] "PC17"                      "PC18"
#[33] "PC19"                      "PC20"
#[35] "PC21"                      "PC22"
## clean expression
geneExprs = log2(recount::getRPKM(rse_gene, "Length")+1)

### regress out after variable 5 (protect 1,2,3,4, 5)
geneExprsClean = cleaningY(geneExprs, modQsva, P=5)

geneExprsClean_sex = cleaningY(geneExprs, modQsva, P=6)

#########################
## get power
powers <- c(1:10, seq(from = 12, to=20, by=2))
sftthresh1 <- pickSoftThreshold(t(geneExprsClean), powerVector = powers,
                               networkType = "signed", verbose = 5)
sftthresh1$powerEstimate
#10
save(sftthresh1, file = "rdas/power_object.rda")

## run wgcna
net = blockwiseModules(t(geneExprsClean), power = sftthresh1$powerEstimate,
                            networkType = "signed", minModuleSize = 30,corType="bicor",
                            reassignThreshold = 0, mergeCutHeight = 0.25,
                            numericLabels = TRUE, pamRespectsDendro = FALSE,
                            saveTOMs = TRUE, verbose = 5, maxBlockSize = 30000,
                            saveTOMFileBase = "rdas/wgcna_signed_TOM")
fNames = rownames(geneExprs)

### saved in later command
#save(net, fNames, file = "rdas/constructed_network_signed_bicor.rda")

########################
## by region remove 
modRegion =  model.matrix(~Dx + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 +
	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr,
	data=colData(rse_gene))

colnames(modRegion)
# [1] "(Intercept)"       "DxControl"         "AgeDeath"
# [4] "SexM"              "snpPC1"            "snpPC2"
# [7] "snpPC3"            "mitoRate"          "rRNA_rate"
# [10] "totalAssignedGene" "RIN"

###jaffe lab function  splits by brain region 
rIndexes = splitit(rse_gene$BrainRegion)

### output of splitit 2 list 
## clean by region ## changed P=3 to P=2 since not clear why to protect AgeDeath
### loop for each brain region


geneExprs_list = mclapply(rIndexes, function(ii) {
##extract the deg exp for each brain region samples
  	degExprs = log2(assays(cov_rse[,ii])$count+1)
#run number of svs for each brain region samples  	
	k = num.sv(degExprs, modRegion[ii,])
#combine model terms/matrix with qSVs *(prcom is a principal component function)l PC on degradation data 	
	m = cbind(modRegion[ii,], prcomp(t(degExprs))$x[,1:k])
#regress out all the effects for gene expression except for intercept and DxControl for each brain region  	
	cleaningY(geneExprs[,ii], m, P=2)
},mc.cores=2)

## threshold
thresh_list = mclapply(geneExprs_list, function(y) {
	pickSoftThreshold(t(y), powerVector = powers,
                               networkType = "signed", verbose = 5)
},mc.cores=2)

## networks
net_list = lapply(1:2, function(i) {
	blockwiseModules(t(geneExprs_list[[i]]), power = thresh_list[[i]]$powerEstimate,
                            networkType = "signed", minModuleSize = 30,corType="bicor",
                            reassignThreshold = 0, mergeCutHeight = 0.25,
                            numericLabels = TRUE, pamRespectsDendro = FALSE,
                            saveTOMs = TRUE, verbose = 5, maxBlockSize = 30000,
                            saveTOMFileBase = paste0("rdas/wgcna_signed_TOM_region",
								names(rIndexes)[i]))
})
save(net_list, net, fNames, file = "rdas/constructed_network_signed_bicor.rda")


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "2020-04-06 14:29:41 EDT"
# 
# 
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 3.6.1 Patched (2019-10-31 r77350)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2020-04-06
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date       lib source
# acepack                1.4.1    2016-10-29 [2] CRAN (R 3.6.1)
# annotate               1.64.0   2019-10-29 [2] Bioconductor
# AnnotationDbi          1.48.0   2019-10-29 [2] Bioconductor
# askpass                1.1      2019-01-13 [2] CRAN (R 3.6.1)
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 3.6.1)
# backports              1.1.5    2019-10-02 [2] CRAN (R 3.6.1)
# base64enc              0.1-3    2015-07-28 [2] CRAN (R 3.6.1)
# Biobase              * 2.46.0   2019-10-29 [2] Bioconductor
# BiocFileCache          1.10.2   2019-11-08 [2] Bioconductor
# BiocGenerics         * 0.32.0   2019-10-29 [2] Bioconductor
# BiocParallel         * 1.20.1   2019-12-21 [2] Bioconductor
# biomaRt                2.42.0   2019-10-29 [2] Bioconductor
# Biostrings             2.54.0   2019-10-29 [2] Bioconductor
# bit                    1.1-15.2 2020-02-10 [2] CRAN (R 3.6.1)
# bit64                  0.9-7    2017-05-08 [2] CRAN (R 3.6.1)
# bitops                 1.0-6    2013-08-17 [2] CRAN (R 3.6.1)
# blob                   1.2.1    2020-01-20 [2] CRAN (R 3.6.1)
# BSgenome               1.54.0   2019-10-29 [2] Bioconductor
# bumphunter             1.28.0   2019-10-29 [2] Bioconductor
# checkmate              2.0.0    2020-02-06 [2] CRAN (R 3.6.1)
# cli                    2.0.2    2020-02-28 [1] CRAN (R 3.6.1)
# cluster                2.1.0    2019-06-19 [3] CRAN (R 3.6.1)
# codetools              0.2-16   2018-12-24 [3] CRAN (R 3.6.1)
# colorout             * 1.2-2    2020-04-02 [1] Github (jalvesaq/colorout@726d681)
# colorspace             1.4-1    2019-03-18 [2] CRAN (R 3.6.1)
# crayon                 1.3.4    2017-09-16 [2] CRAN (R 3.6.1)
# curl                   4.3      2019-12-02 [2] CRAN (R 3.6.1)
# data.table             1.12.8   2019-12-09 [2] CRAN (R 3.6.1)
# DBI                    1.1.0    2019-12-15 [2] CRAN (R 3.6.1)
# dbplyr                 1.4.2    2019-06-17 [2] CRAN (R 3.6.1)
# DelayedArray         * 0.12.2   2020-01-06 [2] Bioconductor
# derfinder              1.20.0   2019-10-29 [2] Bioconductor
# derfinderHelper        1.20.0   2019-10-29 [2] Bioconductor
# digest                 0.6.25   2020-02-23 [1] CRAN (R 3.6.1)
# doParallel             1.0.15   2019-08-02 [2] CRAN (R 3.6.1)
# doRNG                  1.8.2    2020-01-27 [2] CRAN (R 3.6.1)
# downloader             0.4      2015-07-09 [2] CRAN (R 3.6.1)
# dplyr                  0.8.4    2020-01-31 [2] CRAN (R 3.6.1)
# dynamicTreeCut       * 1.63-1   2016-03-11 [1] CRAN (R 3.6.1)
# ellipsis               0.3.0    2019-09-20 [2] CRAN (R 3.6.1)
# fansi                  0.4.1    2020-01-08 [2] CRAN (R 3.6.1)
# fastcluster          * 1.1.25   2018-06-07 [2] CRAN (R 3.6.1)
# foreach                1.4.8    2020-02-09 [2] CRAN (R 3.6.1)
# foreign                0.8-72   2019-08-02 [3] CRAN (R 3.6.1)
# Formula                1.2-3    2018-05-03 [2] CRAN (R 3.6.1)
# genefilter           * 1.68.0   2019-10-29 [2] Bioconductor
# GenomeInfoDb         * 1.22.0   2019-10-29 [2] Bioconductor
# GenomeInfoDbData       1.2.2    2019-10-28 [2] Bioconductor
# GenomicAlignments      1.22.1   2019-11-12 [2] Bioconductor
# GenomicFeatures        1.38.2   2020-02-15 [2] Bioconductor
# GenomicFiles           1.22.0   2019-10-29 [2] Bioconductor
# GenomicRanges        * 1.38.0   2019-10-29 [2] Bioconductor
# GEOquery               2.54.1   2019-11-18 [2] Bioconductor
# ggplot2                3.2.1    2019-08-10 [2] CRAN (R 3.6.1)
# glue                   1.3.2    2020-03-12 [1] CRAN (R 3.6.1)
# GO.db                  3.10.0   2019-10-28 [2] Bioconductor
# googledrive            1.0.0    2019-08-19 [1] CRAN (R 3.6.1)
# gridExtra              2.3      2017-09-09 [2] CRAN (R 3.6.1)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 3.6.1)
# Hmisc                  4.3-1    2020-02-07 [2] CRAN (R 3.6.1)
# hms                    0.5.3    2020-01-08 [2] CRAN (R 3.6.1)
# htmlTable              1.13.3   2019-12-04 [2] CRAN (R 3.6.1)
# htmltools              0.4.0    2019-10-04 [2] CRAN (R 3.6.1)
# htmlwidgets            1.5.1    2019-10-08 [2] CRAN (R 3.6.1)
# httr                   1.4.1    2019-08-05 [2] CRAN (R 3.6.1)
# impute                 1.60.0   2019-10-29 [2] Bioconductor
# IRanges              * 2.20.2   2020-01-13 [2] Bioconductor
# iterators              1.0.12   2019-07-26 [2] CRAN (R 3.6.1)
# jaffelab             * 0.99.30  2020-04-02 [1] Github (LieberInstitute/jaffelab@42637ff)
# jpeg                   0.1-8.1  2019-10-24 [2] CRAN (R 3.6.1)
# jsonlite               1.6.1    2020-02-02 [2] CRAN (R 3.6.1)
# knitr                  1.28     2020-02-06 [2] CRAN (R 3.6.1)
# lattice                0.20-38  2018-11-04 [3] CRAN (R 3.6.1)
# latticeExtra           0.6-29   2019-12-19 [2] CRAN (R 3.6.1)
# lazyeval               0.2.2    2019-03-15 [2] CRAN (R 3.6.1)
# lifecycle              0.2.0    2020-03-06 [1] CRAN (R 3.6.1)
# limma                  3.42.2   2020-02-03 [2] Bioconductor
# locfit                 1.5-9.1  2013-04-20 [2] CRAN (R 3.6.1)
# magrittr               1.5      2014-11-22 [2] CRAN (R 3.6.1)
# Matrix                 1.2-17   2019-03-22 [3] CRAN (R 3.6.1)
# matrixStats          * 0.55.0   2019-09-07 [2] CRAN (R 3.6.1)
# memoise                1.1.0    2017-04-21 [2] CRAN (R 3.6.1)
# mgcv                 * 1.8-30   2019-10-24 [3] CRAN (R 3.6.1)
# munsell                0.5.0    2018-06-12 [2] CRAN (R 3.6.1)
# nlme                 * 3.1-141  2019-08-01 [3] CRAN (R 3.6.1)
# nnet                   7.3-12   2016-02-02 [3] CRAN (R 3.6.1)
# openssl                1.4.1    2019-07-18 [2] CRAN (R 3.6.1)
# pillar                 1.4.3    2019-12-20 [2] CRAN (R 3.6.1)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 3.6.1)
# plyr                   1.8.5    2019-12-10 [2] CRAN (R 3.6.1)
# png                    0.1-7    2013-12-03 [2] CRAN (R 3.6.1)
# preprocessCore         1.48.0   2019-10-29 [2] Bioconductor
# prettyunits            1.1.1    2020-01-24 [2] CRAN (R 3.6.1)
# progress               1.2.2    2019-05-16 [2] CRAN (R 3.6.1)
# purrr                  0.3.3    2019-10-18 [2] CRAN (R 3.6.1)
# qvalue                 2.18.0   2019-10-29 [2] Bioconductor
# R6                     2.4.1    2019-11-12 [2] CRAN (R 3.6.1)
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 3.6.1)
# rappdirs               0.3.1    2016-03-28 [2] CRAN (R 3.6.1)
# RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 3.6.1)
# Rcpp                   1.0.3    2019-11-08 [2] CRAN (R 3.6.1)
# RCurl                  1.98-1.1 2020-01-19 [2] CRAN (R 3.6.1)
# readr                  1.3.1    2018-12-21 [2] CRAN (R 3.6.1)
# recount              * 1.12.1   2019-11-06 [2] Bioconductor
# rentrez                1.2.2    2019-05-02 [2] CRAN (R 3.6.1)
# reshape2               1.4.3    2017-12-11 [2] CRAN (R 3.6.1)
# rlang                  0.4.5    2020-03-01 [1] CRAN (R 3.6.1)
# rngtools               1.5      2020-01-23 [2] CRAN (R 3.6.1)
# rpart                  4.1-15   2019-04-12 [3] CRAN (R 3.6.1)
# Rsamtools              2.2.2    2020-02-11 [2] Bioconductor
# RSQLite                2.2.0    2020-01-07 [2] CRAN (R 3.6.1)
# rstudioapi             0.11     2020-02-07 [2] CRAN (R 3.6.1)
# rtracklayer            1.46.0   2019-10-29 [2] Bioconductor
# S4Vectors            * 0.24.3   2020-01-18 [2] Bioconductor
# scales                 1.1.0    2019-11-18 [2] CRAN (R 3.6.1)
# segmented              1.1-0    2019-12-10 [2] CRAN (R 3.6.1)
# sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 3.6.1)
# stringi                1.4.6    2020-02-17 [2] CRAN (R 3.6.1)
# stringr                1.4.0    2019-02-10 [2] CRAN (R 3.6.1)
# SummarizedExperiment * 1.16.1   2019-12-19 [2] Bioconductor
# survival               3.1-8    2019-12-03 [2] CRAN (R 3.6.1)
# sva                  * 3.34.0   2019-10-29 [2] Bioconductor
# tibble                 3.0.0    2020-03-30 [1] CRAN (R 3.6.1)
# tidyr                  1.0.2    2020-01-24 [2] CRAN (R 3.6.1)
# tidyselect             1.0.0    2020-01-27 [2] CRAN (R 3.6.1)
# VariantAnnotation      1.32.0   2019-10-29 [2] Bioconductor
# vctrs                  0.2.4    2020-03-10 [1] CRAN (R 3.6.1)
# WGCNA                * 1.69     2020-02-28 [1] CRAN (R 3.6.1)
# withr                  2.1.2    2018-03-15 [2] CRAN (R 3.6.1)
# xfun                   0.12     2020-01-13 [2] CRAN (R 3.6.1)
# XML                    3.99-0.3 2020-01-20 [2] CRAN (R 3.6.1)
# xml2                   1.2.2    2019-08-09 [2] CRAN (R 3.6.1)
# xtable                 1.8-4    2019-04-21 [2] CRAN (R 3.6.1)
# XVector                0.26.0   2019-10-29 [2] Bioconductor
# zlibbioc               1.32.0   2019-10-29 [2] Bioconductor
# 
# [1] /users/fgoes/R/3.6.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/library
# 


