## Load plink before starting R
# module load plink/1.90b6.6
# R

## Now run R code
library("data.table")
library("SummarizedExperiment")
library("here")
library("recount")
library("sva")
library("sessioninfo")
library("tidyverse")
library("getopt")

spec <- matrix(
    c('region', 'r', 1, 'character', 'Either Amygdala or SACC'),
    byrow = TRUE,
    ncol = 5
)
opt <- getopt(spec)

# opt$region <- "Amygdala"

if (tolower(opt$region) == "sacc") {
    opt$region = "sACC"
} else if (tolower(opt$region) == "amygdala" ||
           tolower(opt$region) == "amyg") {
    opt$region = "Amygdala"
} else {
    cat(getopt(spec, usage = TRUE))
    q(status = 1)
}

## Show the options used
message(paste(Sys.time(), "options used"))
print(opt)

dir.create(paste0(opt$region, "_rda"), showWarnings = FALSE)

## To avoid issues with running this code on qsub
data.table::setDTthreads(threads = 1)

## Find the samples for this project
load(here::here("exprs_cutoff", "rse_gene.Rdata"))

## 99 brains only have 1 sample either in one region or another
## 540 amyg samples
## 551 sACC samples
## 588 MDD
## 503 BP

# subset rse_gene to either Amygdala or sACC, whichever the user selected
rse_gene <- rse_gene[, colData(rse_gene)$BrainRegion == opt$region]

stopifnot(length(unique(rse_gene$BrNum)) == ncol(rse_gene))

# Load snpPCs
load(here::here("genotype_data", "goesHyde_bipolarMdd_Genotypes_mds.rda"))

## For converting BrNum's into numbers
brnumerical <- function(x) {
    as.integer(gsub("Br|_.*", "", x))
}

# Read in original PLINK files
libd_bfile <-
    here::here("genotype_data", "mdd_bpd", "maf01", "mdd_bpd_maf01.rsid")

# Cross referencing to retrieve BrNums
cross_ref <-
    fread(
        "/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/merged_batches_topmed/usable_genotypes/maf01/Genotype2422.csv",
        skip = 1,
        drop = 1
    )
colnames(cross_ref) <- c("PLINK_ID", "BrNum")

new_plink_loc <-
    here::here("twas_both", "filter_data", "reformat_genotype")

dir.create(new_plink_loc, showWarnings = FALSE)

# copying the genotype data
if (!file.exists(here(
    "twas_both",
    "filter_data",
    "reformat_genotype",
    "mdd_bpd_maf01.rsid.bed"
))) {
    system(paste0("cp ", libd_bfile, "* ", new_plink_loc, "/"))
}

## Read the LIBD fam data
libd_fam <- fread(
    paste0(libd_bfile, ".fam"),
    col.names = c(
        "famid",
        "w_famid",
        "w_famid_fa",
        "w_famid_mo",
        "sex_code",
        "phenotype"
    )
)

libd_fam$plink_key <- paste0(libd_fam$famid, "_", libd_fam$w_famid)

libd_fam$BrNum <-
    cross_ref[match(libd_fam$plink_key, cross_ref$PLINK_ID), ]$BrNum

libd_fam_save <- libd_fam %>%
    mutate(
        famid = BrNum,
        w_famid = BrNum,
        plink_key = NULL,
        BrNum = NULL
    )

libd_bfile <-
    here::here("twas_both",
               "filter_data",
               "reformat_genotype",
               "mdd_bpd_maf01.rsid")

fwrite(
    libd_fam_save,
    file = paste0(libd_bfile, ".fam"),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
)

libd_fam$brnumerical <- brnumerical(libd_fam$BrNum)
setkey(libd_fam, "brnumerical")

## Filter the LIBD data to the one specific to this project
message(paste(Sys.time(), "processing", opt$region))
samp_file <- paste0("samples_to_extract_", opt$region, ".txt")

## Which samples have genotype data and MDS data?
samples_in_all <- intersect(intersect(brnumerical(rse_gene$BrNum), libd_fam$brnumerical),
                            brnumerical(rownames(mds)))

################################
## Subset and save all key files
################################
rse_gene <-
    rse_gene[, brnumerical(rse_gene$BrNum) %in% samples_in_all]

## Match rse_gene to mds
m_to_mds <-
    match(brnumerical(rse_gene$BrNum), brnumerical(rownames(mds)))
stopifnot(all(!is.na(m_to_mds)))
mds <- mds[m_to_mds,]

## Append mds to colData
colData(rse_gene) <- cbind(colData(rse_gene), mds)

## Compute RPKM
assays(rse_gene)$RPKM <- getRPKM(rse_gene, "Length")

## Compute gene PCs
message(Sys.time(), " computing gene PCs on log2(RPKM + 1)")
pcaGene <- prcomp(t(log2(assays(rse_gene)$RPKM + 1)))
save(pcaGene, file = paste0(opt$region, "_rda/pcaGene.Rdata"))

message(Sys.time(), " determine how many gene PCs to adjust for")
mod <-
    model.matrix( ~ Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, data = colData(rse_gene))
kGene <- num.sv(log2(assays(rse_gene)$RPKM + 1), mod)
stopifnot(kGene > 0)
genePCs <- pcaGene$x[, seq_len(kGene)]
save(genePCs, file = paste0(opt$region, "_rda/genePCs.Rdata"))

## Add gene PCs to rse_gene
colData(rse_gene) <- cbind(colData(rse_gene), genePCs)

## Save for later
save(
    rse_gene,
    file = paste0(
        opt$region,
        "_rda/",
        opt$region,
        "_hg38_rseGene_rawCounts_allSamples_n",
        ncol(rse_gene),
        ".Rdata"
    )
)

## Now extract the genotype data too
filter_m <- match(brnumerical(rse_gene$BrNum), libd_fam$brnumerical)
stopifnot(all(!is.na(filter_m)))
fwrite(libd_fam[filter_m, 1:2],
       ## can be more involved
       file = samp_file,
       sep = "\t",
       col.names = FALSE)
newbfile_root <- paste0("mdd_bpd_maf01.rsid", opt$region)

dir.create(paste0(opt$region, "_duplicate_snps_bim"), showWarnings = FALSE)
newbfile <-
    here::here(
        "twas_both",
        "filter_data",
        paste0(opt$region, "_duplicate_snps_bim"),
        paste0(newbfile_root,
               "_duplicateSNPs")
    )

## Extract
message(paste(Sys.time(), "running bfile extract for", newbfile))
system(
    paste(
        "plink --bfile",
        libd_bfile,
        "--keep",
        samp_file,
        "--make-bed --out",
        newbfile,
        " --memory 100000 --biallelic-only"
    )
)

## Check that we have the right data
newbfile_fam <- fread(
    paste0(newbfile, ".fam"),
    col.names = c(
        "famid",
        "w_famid",
        "w_famid_fa",
        "w_famid_mo",
        "sex_code",
        "phenotype"
    )
)
check_m <-
    match(brnumerical(newbfile_fam$famid),
          brnumerical(colData(rse_gene)$BrNum))
stopifnot(all(!is.na(check_m)))


## Re-run but now make the SNV names unique
dir.create(paste0(opt$region, "_unique_snps_bim"), showWarnings = FALSE)
newbfile_unique <-
    here::here(
        "twas_both",
        "filter_data",
        paste0(opt$region, "_unique_snps_bim"),
        paste0(newbfile_root,
               "_uniqueSNPs")
    )

## Extract again (could also copy and rename, but it's fast this way)
message(paste(Sys.time(), "running bfile extract for", newbfile_unique))
system(
    paste(
        "plink --bfile",
        libd_bfile,
        "--keep",
        samp_file,
        "--make-bed --out",
        newbfile_unique,
        " --memory 100000 --biallelic-only"
    )
)


message(paste(Sys.time(), "reading the bim file", newbfile_unique))
bim <- fread(
    paste0(newbfile_unique, ".bim"),
    col.names = c("chr", "snp", "position", "basepair", "allele1", "allele2")
)

table(duplicated(bim$snp))
#    FALSE     TRUE
# 10943065    44114

## Make names unique
message(Sys.time(), " making the variant names unique")
bim$snp <- make.names(bim$snp, unique = TRUE)
stopifnot(all(!duplicated(bim$snp)))

## Ovewrite the PLINK bim file
fwrite(
    bim,
    file = paste0(newbfile_unique, ".bim"),
    sep = " ",
    col.names = FALSE
)


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
# date     2021-09-22                                 
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date       lib source        
# annotate               1.70.0   2021-05-19 [2] Bioconductor  
# AnnotationDbi          1.54.1   2021-06-08 [2] Bioconductor  
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# backports              1.2.1    2020-12-09 [2] CRAN (R 4.1.0)
# base64enc              0.1-3    2015-07-28 [2] CRAN (R 4.1.0)
# Biobase              * 2.52.0   2021-05-19 [2] Bioconductor  
# BiocFileCache          2.0.0    2021-05-19 [2] Bioconductor  
# BiocGenerics         * 0.38.0   2021-05-19 [2] Bioconductor  
# BiocIO                 1.2.0    2021-05-19 [2] Bioconductor  
# BiocParallel         * 1.26.1   2021-07-04 [2] Bioconductor  
# biomaRt                2.48.2   2021-07-01 [2] Bioconductor  
# Biostrings             2.60.2   2021-08-05 [2] Bioconductor  
# bit                    4.0.4    2020-08-04 [2] CRAN (R 4.1.0)
# bit64                  4.0.5    2020-08-30 [2] CRAN (R 4.1.0)
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# blob                   1.2.2    2021-07-23 [2] CRAN (R 4.1.0)
# broom                  0.7.9    2021-07-27 [2] CRAN (R 4.1.0)
# BSgenome               1.60.0   2021-05-19 [2] Bioconductor  
# bumphunter             1.34.0   2021-05-19 [2] Bioconductor  
# cachem                 1.0.5    2021-05-15 [2] CRAN (R 4.1.0)
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.1.0)
# checkmate              2.0.0    2020-02-06 [2] CRAN (R 4.1.0)
# cli                    3.0.1    2021-07-17 [2] CRAN (R 4.1.0)
# cluster                2.1.2    2021-04-17 [3] CRAN (R 4.1.0)
# codetools              0.2-18   2020-11-04 [2] CRAN (R 4.1.0)
# colorspace             2.0-2    2021-06-24 [2] CRAN (R 4.1.0)
# crayon                 1.4.1    2021-02-08 [2] CRAN (R 4.1.0)
# curl                   4.3.2    2021-06-23 [2] CRAN (R 4.1.0)
# data.table           * 1.14.0   2021-02-21 [2] CRAN (R 4.1.0)
# DBI                    1.1.1    2021-01-15 [2] CRAN (R 4.1.0)
# dbplyr                 2.1.1    2021-04-06 [2] CRAN (R 4.1.0)
# DelayedArray           0.18.0   2021-05-19 [2] Bioconductor  
# derfinder              1.26.0   2021-05-19 [2] Bioconductor  
# derfinderHelper        1.26.0   2021-05-19 [2] Bioconductor  
# digest                 0.6.27   2020-10-24 [2] CRAN (R 4.1.0)
# doRNG                  1.8.2    2020-01-27 [2] CRAN (R 4.1.0)
# downloader             0.4      2015-07-09 [2] CRAN (R 4.1.0)
# dplyr                * 1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
# edgeR                  3.34.0   2021-05-19 [2] Bioconductor  
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  0.5.0    2021-05-25 [2] CRAN (R 4.1.0)
# fastmap                1.1.0    2021-01-25 [2] CRAN (R 4.1.0)
# filelock               1.0.2    2018-10-05 [2] CRAN (R 4.1.0)
# forcats              * 0.5.1    2021-01-27 [2] CRAN (R 4.1.0)
# foreach                1.5.1    2020-10-15 [2] CRAN (R 4.1.0)
# foreign                0.8-81   2020-12-22 [3] CRAN (R 4.1.0)
# Formula                1.2-4    2020-10-16 [2] CRAN (R 4.1.0)
# fs                     1.5.0    2020-07-31 [2] CRAN (R 4.1.0)
# genefilter           * 1.74.0   2021-05-19 [2] Bioconductor  
# generics               0.1.0    2020-10-31 [2] CRAN (R 4.1.0)
# GenomeInfoDb         * 1.28.1   2021-07-01 [2] Bioconductor  
# GenomeInfoDbData       1.2.6    2021-05-11 [2] Bioconductor  
# GenomicAlignments      1.28.0   2021-05-19 [2] Bioconductor  
# GenomicFeatures        1.44.0   2021-05-19 [2] Bioconductor  
# GenomicFiles           1.28.0   2021-05-19 [2] Bioconductor  
# GenomicRanges        * 1.44.0   2021-05-19 [2] Bioconductor  
# GEOquery               2.60.0   2021-05-19 [2] Bioconductor  
# getopt               * 1.20.3   2019-03-22 [2] CRAN (R 4.1.0)
# ggplot2              * 3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
# glue                   1.4.2    2020-08-27 [2] CRAN (R 4.1.0)
# gridExtra              2.3      2017-09-09 [2] CRAN (R 4.1.0)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# haven                  2.4.3    2021-08-04 [2] CRAN (R 4.1.0)
# here                 * 1.0.1    2020-12-13 [1] CRAN (R 4.1.0)
# Hmisc                  4.5-0    2021-02-28 [2] CRAN (R 4.1.0)
# hms                    1.1.0    2021-05-17 [2] CRAN (R 4.1.0)
# htmlTable              2.2.1    2021-05-18 [2] CRAN (R 4.1.0)
# htmltools              0.5.1.1  2021-01-22 [2] CRAN (R 4.1.0)
# htmlwidgets            1.5.3    2020-12-10 [2] CRAN (R 4.1.0)
# httr                   1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
# IRanges              * 2.26.0   2021-05-19 [2] Bioconductor  
# iterators              1.0.13   2020-10-15 [2] CRAN (R 4.1.0)
# jpeg                   0.1-9    2021-07-24 [2] CRAN (R 4.1.0)
# jsonlite               1.7.2    2020-12-09 [2] CRAN (R 4.1.0)
# KEGGREST               1.32.0   2021-05-19 [2] Bioconductor  
# knitr                  1.33     2021-04-24 [2] CRAN (R 4.1.0)
# lattice                0.20-44  2021-05-02 [3] CRAN (R 4.1.0)
# latticeExtra           0.6-29   2019-12-19 [2] CRAN (R 4.1.0)
# lifecycle              1.0.0    2021-02-15 [2] CRAN (R 4.1.0)
# limma                  3.48.3   2021-08-10 [2] Bioconductor  
# locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
# lubridate              1.7.10   2021-02-26 [2] CRAN (R 4.1.0)
# magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
# Matrix                 1.3-4    2021-06-01 [3] CRAN (R 4.1.0)
# MatrixGenerics       * 1.4.2    2021-08-08 [2] Bioconductor  
# matrixStats          * 0.60.0   2021-07-26 [2] CRAN (R 4.1.0)
# memoise                2.0.0    2021-01-26 [2] CRAN (R 4.1.0)
# mgcv                 * 1.8-36   2021-06-01 [3] CRAN (R 4.1.0)
# modelr                 0.1.8    2020-05-19 [2] CRAN (R 4.1.0)
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# nlme                 * 3.1-152  2021-02-04 [3] CRAN (R 4.1.0)
# nnet                   7.3-16   2021-05-03 [3] CRAN (R 4.1.0)
# pillar                 1.6.2    2021-07-29 [2] CRAN (R 4.1.0)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# plyr                   1.8.6    2020-03-03 [2] CRAN (R 4.1.0)
# png                    0.1-7    2013-12-03 [2] CRAN (R 4.1.0)
# prettyunits            1.1.1    2020-01-24 [2] CRAN (R 4.1.0)
# progress               1.2.2    2019-05-16 [2] CRAN (R 4.1.0)
# purrr                * 0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# qvalue                 2.24.0   2021-05-19 [2] Bioconductor  
# R6                     2.5.0    2020-10-28 [2] CRAN (R 4.1.0)
# rappdirs               0.3.3    2021-01-31 [2] CRAN (R 4.1.0)
# RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
# Rcpp                   1.0.7    2021-07-07 [2] CRAN (R 4.1.0)
# RCurl                  1.98-1.3 2021-03-16 [2] CRAN (R 4.1.0)
# readr                * 2.0.1    2021-08-10 [2] CRAN (R 4.1.0)
# readxl                 1.3.1    2019-03-13 [2] CRAN (R 4.1.0)
# recount              * 1.18.1   2021-08-10 [2] Bioconductor  
# rentrez                1.2.3    2020-11-10 [2] CRAN (R 4.1.0)
# reprex                 2.0.1    2021-08-05 [2] CRAN (R 4.1.0)
# reshape2               1.4.4    2020-04-09 [1] CRAN (R 4.1.0)
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
# rvest                  1.0.1    2021-07-26 [2] CRAN (R 4.1.0)
# S4Vectors            * 0.30.0   2021-05-19 [2] Bioconductor  
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
# sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.1.0)
# stringi                1.7.3    2021-07-16 [2] CRAN (R 4.1.0)
# stringr              * 1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment * 1.22.0   2021-05-19 [2] Bioconductor  
# survival               3.2-12   2021-08-13 [3] CRAN (R 4.1.0)
# sva                  * 3.40.0   2021-05-19 [2] Bioconductor  
# tibble               * 3.1.3    2021-07-23 [2] CRAN (R 4.1.0)
# tidyr                * 1.1.3    2021-03-03 [2] CRAN (R 4.1.0)
# tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
# tidyverse            * 1.3.1    2021-04-15 [2] CRAN (R 4.1.0)
# tzdb                   0.1.2    2021-07-20 [2] CRAN (R 4.1.0)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# VariantAnnotation      1.38.0   2021-05-19 [2] Bioconductor  
# vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
# withr                  2.4.2    2021-04-18 [2] CRAN (R 4.1.0)
# xfun                   0.25     2021-08-06 [2] CRAN (R 4.1.0)
# XML                    3.99-0.6 2021-03-16 [2] CRAN (R 4.1.0)
# xml2                   1.3.2    2020-04-23 [2] CRAN (R 4.1.0)
# xtable                 1.8-4    2019-04-21 [2] CRAN (R 4.1.0)
# XVector                0.32.0   2021-05-19 [2] Bioconductor  
# yaml                   2.2.1    2020-02-01 [2] CRAN (R 4.1.0)
# zlibbioc               1.38.0   2021-05-19 [2] Bioconductor  
# 
# [1] /users/aseyedia/R/4.1
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/library

# PLINK v1.90b6.6 64-bit (10 Oct 2018)

system("plink --version")

