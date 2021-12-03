
library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(sessioninfo)
library(here)

##### Load rse data, examine ####
#load objects
load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)
pd = colData(rse_gene)

table(pd$Experiment)
# psychENCODE_MDD  psychENCODE_BP 
# 588             503 

## load degradation data
load(here("data","degradation_rse_MDDseq_BiPSeq_BothRegions.Rdata"), verbose = TRUE)

##### get qSVs ####
message("Basic Model")
modJoint = model.matrix(~PrimaryDx*BrainRegion + AgeDeath + Sex + 
                          snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
                          mitoRate + rRNA_rate + totalAssignedGene + RIN + abs(ERCCsumLogErr),
                        data=colData(rse_gene))
colnames(modJoint)

# [1] "(Intercept)"                      "PrimaryDxControl"                 "PrimaryDxBipolar"                
# [4] "BrainRegionsACC"                  "AgeDeath"                         "SexM"                            
# [7] "snpPC1"                           "snpPC2"                           "snpPC3"                          
# [10] "snpPC4"                           "snpPC5"                           "mitoRate"                        
# [13] "rRNA_rate"                        "totalAssignedGene"                "RIN"                             
# [16] "abs(ERCCsumLogErr)"               "PrimaryDxControl:BrainRegionsACC" "PrimaryDxBipolar:BrainRegionsACC"
message("\nBasic Model + Cell Fractions")
modJoint_deconvo = model.matrix(~PrimaryDx*BrainRegion + AgeDeath + Sex + 
                          snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
                          mitoRate + rRNA_rate + totalAssignedGene + RIN + abs(ERCCsumLogErr) +
                          Astro + Endo + Macro + Micro + Mural + Oligo + OPC + Tcell + Excit,
                        data=colData(rse_gene))

colnames(modJoint_deconvo)
# [1] "(Intercept)"                      "PrimaryDxControl"                 "PrimaryDxBipolar"                
# [4] "BrainRegionsACC"                  "AgeDeath"                         "SexM"                            
# [7] "snpPC1"                           "snpPC2"                           "snpPC3"                          
# [10] "snpPC4"                           "snpPC5"                           "mitoRate"                        
# [13] "rRNA_rate"                        "totalAssignedGene"                "RIN"                             
# [16] "ERCCsumLogErr"                    "Astro"                            "Endo"                            
# [19] "Macro"                            "Micro"                            "Mural"                           
# [22] "Oligo"                            "OPC"                              "Tcell"                           
# [25] "Excit"                            "PrimaryDxControl:BrainRegionsACC" "PrimaryDxBipolar:BrainRegionsACC"

message("\nModel with Experiment Term")
modExp = model.matrix(~PrimaryDx*BrainRegion + AgeDeath + Sex + Experiment +
                        snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
                        mitoRate + rRNA_rate + totalAssignedGene + RIN + abs(ERCCsumLogErr),
                      data=colData(rse_gene))
colnames(modExp)

## Save Models
save(modJoint, modJoint_deconvo, file = here("differential_expression","data","differental_modJoint.Rdata"))

degExprs = log2(assays(cov_rse)$count+1)
k = num.sv(degExprs, modJoint)
message("k=", k)
# k=26

k_deconvo = num.sv(degExprs, modJoint_deconvo)
message("k deconvo = ", k_deconvo)
# k deconvo = 26

k_exp = num.sv(degExprs, modExp)
message("k Experiment = ", k_exp)
# k Experiment = 25

qSV_mat = prcomp(t(degExprs))$x[,1:k]

## Save qSVmat
save(qSV_mat, file = here("differential_expression","data","qSV_mat.Rdata"))

## Explore variance explained
varExplQsva = jaffelab::getPcaVars(prcomp(t(degExprs)))
varExplQsva[1:k]

# [1] 67.700  4.860  3.010  1.990  1.430  1.210  1.120  0.803  0.778  0.677
# [11]  0.574  0.521  0.487  0.444  0.347  0.342  0.329  0.300  0.246  0.237
# [21]  0.230  0.223  0.205  0.196  0.190  0.174

sum(varExplQsva[1:k]) 
# [1] 88.623

#sgejobs::job_single('qSV_calculation', create_shell = TRUE, memory = '80G', command = "Rscript qSV_calculations.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# [1] "2021-12-03 11:10:25 EST"
# user  system elapsed 
# 251.827   1.404 254.231 
# ─ Session info  ──────────────────────────────────────────────────────────────────────────────────────────────────────
# hash: person in suit levitating: dark skin tone, notebook with decorative cover, briefcase
# 
# setting  value
# version  R version 4.1.2 Patched (2021-11-04 r81138)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2021-12-03
# pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date (UTC) lib source
# annotate               1.72.0   2021-10-26 [2] Bioconductor
# AnnotationDbi          1.56.2   2021-11-09 [2] Bioconductor
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# Biobase              * 2.54.0   2021-10-26 [2] Bioconductor
# BiocGenerics         * 0.40.0   2021-10-26 [2] Bioconductor
# BiocParallel         * 1.28.2   2021-11-25 [2] Bioconductor
# Biostrings             2.62.0   2021-10-26 [2] Bioconductor
# bit                    4.0.4    2020-08-04 [2] CRAN (R 4.1.0)
# bit64                  4.0.5    2020-08-30 [2] CRAN (R 4.1.0)
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# blob                   1.2.2    2021-07-23 [2] CRAN (R 4.1.0)
# cachem                 1.0.6    2021-08-19 [2] CRAN (R 4.1.2)
# cli                    3.1.0    2021-10-27 [2] CRAN (R 4.1.2)
# crayon                 1.4.2    2021-10-29 [2] CRAN (R 4.1.2)
# DBI                    1.1.1    2021-01-15 [2] CRAN (R 4.1.0)
# DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
# dplyr                  1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
# edgeR                  3.36.0   2021-10-26 [2] Bioconductor
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  0.5.0    2021-05-25 [2] CRAN (R 4.1.0)
# fastmap                1.1.0    2021-01-25 [2] CRAN (R 4.1.0)
# fs                     1.5.1    2021-11-30 [2] CRAN (R 4.1.2)
# gargle                 1.2.0    2021-07-02 [2] CRAN (R 4.1.0)
# genefilter           * 1.76.0   2021-10-26 [2] Bioconductor
# generics               0.1.1    2021-10-25 [2] CRAN (R 4.1.2)
# GenomeInfoDb         * 1.30.0   2021-10-26 [2] Bioconductor
# GenomeInfoDbData       1.2.7    2021-11-01 [2] Bioconductor
# GenomicRanges        * 1.46.1   2021-11-18 [2] Bioconductor
# glue                   1.5.1    2021-11-30 [2] CRAN (R 4.1.2)
# googledrive            2.0.0    2021-07-08 [2] CRAN (R 4.1.0)
# here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
# httr                   1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
# IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
# jaffelab             * 0.99.31  2021-10-29 [1] Github (LieberInstitute/jaffelab@2cbd55a)
# KEGGREST               1.34.0   2021-10-26 [2] Bioconductor
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lifecycle              1.0.1    2021-09-24 [1] CRAN (R 4.1.1)
# limma                  3.50.0   2021-10-26 [2] Bioconductor
# locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
# magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
# Matrix                 1.3-4    2021-06-01 [3] CRAN (R 4.1.2)
# MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
# matrixStats          * 0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
# memoise                2.0.1    2021-11-26 [2] CRAN (R 4.1.2)
# mgcv                 * 1.8-38   2021-10-06 [3] CRAN (R 4.1.2)
# nlme                 * 3.1-153  2021-09-07 [3] CRAN (R 4.1.2)
# pillar                 1.6.4    2021-10-18 [1] CRAN (R 4.1.1)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# png                    0.1-7    2013-12-03 [2] CRAN (R 4.1.0)
# purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R6                     2.5.1    2021-08-19 [1] CRAN (R 4.1.1)
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.1.1)
# RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
# Rcpp                   1.0.7    2021-07-07 [2] CRAN (R 4.1.0)
# RCurl                  1.98-1.5 2021-09-17 [2] CRAN (R 4.1.2)
# rlang                  0.4.12   2021-10-18 [2] CRAN (R 4.1.2)
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
# RSQLite                2.2.8    2021-08-21 [2] CRAN (R 4.1.2)
# S4Vectors            * 0.32.3   2021-11-21 [2] Bioconductor
# segmented              1.3-4    2021-04-22 [1] CRAN (R 4.1.1)
# sessioninfo          * 1.2.1    2021-11-02 [2] CRAN (R 4.1.2)
# SummarizedExperiment * 1.24.0   2021-10-26 [2] Bioconductor
# survival               3.2-13   2021-08-24 [3] CRAN (R 4.1.2)
# sva                  * 3.42.0   2021-10-26 [2] Bioconductor
# tibble                 3.1.6    2021-11-07 [2] CRAN (R 4.1.2)
# tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
# XML                    3.99-0.8 2021-09-17 [2] CRAN (R 4.1.2)
# xtable                 1.8-4    2019-04-21 [2] CRAN (R 4.1.0)
# XVector                0.34.0   2021-10-26 [2] Bioconductor
# zlibbioc               1.40.0   2021-10-26 [2] Bioconductor
# 
# [1] /users/lhuuki/R/4.1.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# **** Job ends ****
#   Fri Dec  3 11:10:26 EST 2021
