library("readxl")
library("sessioninfo")
library("SummarizedExperiment")
library("jaffelab")
library("dplyr")
library("reshape2")
library("purrr")
library("here")

## style this script
# styler::style_file("build_metadata.R", transformers = styler::tidyverse_style(indent_by = 4))

source("build_metadata_functions.R")

#### load data ####
lims <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/shiny070220.csv")
pd <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/pd_swap.csv", as.is = TRUE)

# add brain weight to pd
brain_weight <- lims %>% select(BrNum, `Brain.Weight..gram.`)

pd_all <- read.csv(here("data","raw_GoesZandi_pd.csv")) %>%
    filter(!overlap | Experiment == "psychENCODE_MDD") # remove overlaps
pd_mdd <- pd_all %>%
    filter(Experiment == "psychENCODE_MDD") %>%
    left_join(brain_weight, by = "BrNum") %>%
    mutate(BrodmannArea = ifelse(BrainRegion == "anterior cingulate cortex", 25, NA))

# add imputed Sex to lims
brain_sentrix <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/brain_sentrix_swap.csv") %>%
    filter(ID %in% pd_mdd$genoSample)

lims <- lims %>% left_join(brain_sentrix %>% select(BrNum, genoSex))

# replace values
gia <- read.csv("genotypeInferredAncestry.csv")

lims <- lims %>%
    left_join(gia) %>%
    mutate(genotypeInferredAncestry = ifelse(is.na(genotypeInferredAncestry), Race, genotypeInferredAncestry)) %>%
    select(-PrimaryDx) %>%
    inner_join(pd_mdd %>% select(PrimaryDx, BrNum) %>% unique(), by = "BrNum")


# fastq info
fastq_info_fn <- "fastq_info_mdd.csv"
if(!file.exists("fastq_info_mdd.csv")){
    fastq_info <- map(mdd_manifest$path, get_fastq_info)
    fastq_info_df <- fastq_info %>%
        as.data.frame() %>%
        t()
    write.csv(fastq_info_df, file = fastq_info_fn)
}else{
    fastq_info_df <-read.csv(fastq_info_fn)
}

fastq_info_df <- cbind(fastq_info_df, mdd_manifest[, "RNum", drop = FALSE])
rownames(fastq_info_df) <- NULL

flow_cell <- fastq_info_df %>%
    select(flow_cell, RNum) %>%
    unique()
pd <- pd %>%
    left_join(flow_cell, by = "RNum") %>%
    filter(Dataset == "psychENCODE_MDD")

# list samples
BrNum_mdd <- unique(pd_mdd$BrNum)
RNum_mdd <- pd_mdd$RNum

## Build genodata table
pd_geno <- pd_all %>%
    select(BrNum, genoSample)%>%
    unique %>%
    left_join(brain_sentrix %>% select(genoSample = ID, BrNum, Batch))

geno_files <- scan(here("genotype_data","manifest.txt"), what="character", sep="\n")
geno_assay <- tibble(file = geno_files,
                     fileFormat = gsub("^.*\\.", "", geno_files)) %>%
    filter(fileFormat == "fam") %>%
    mutate(chip = gsub("_md\\S+.fam","",file))

fam_tables <- map(geno_assay$file, ~read.delim(here("genotype_data",.x), header = FALSE,sep = " "))
fam_df <- bind_rows(fam_tables) %>%
    transmute(genoSample = paste0(V1, "_", V2))
fam_df$file <- rep(geno_assay$file, map_int(fam_tables, nrow))
fam_df <- fam_df %>% left_join(geno_assay)
dim(fam_df)
# [1] 1204    4
fam_df_pd <- fam_df %>% left_join(pd_geno) %>%
    mutate(Batch = ifelse(chip == "topmed", NA, Batch))
# More info here: dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/
fam_df_pd %>% count(chip)
#                        chip   n
# 1               Illumina_1M 125
# 2 Illumina_HumanHap650Yv3_A  25
# 3   Illumina_Omni2.5-8_v1.1 143
# 4   Illumina_Omni2.5-8_v1.3 212
# 5   Illumina_Omni2.5-8_v1.4  23
# 6   Illumina_Omni2.5-8_v1.5  19
# 7       Illumina_Omni5-Quad  55
# 8                    topmed 602
fam_df_pd %>% count(chip == "topmed")
# chip == "topmed"   n
# 1            FALSE 602
# 2             TRUE 602
fam_df_mdd <- fam_df_pd %>% filter(BrNum %in% pd_mdd$BrNum)
## Save annotation file
annotation_fn = "psychENCODE_MDD_genotype_file_annotation.tsv"
fam_df_mdd  %>%
    select(Individual = BrNum, specimenID = genoSample, file) %>%
    write.table(file = annotation_fn, row.names = FALSE, sep = "\t")

genoSample_mdd <- unique(fam_df_mdd$genoSample)

#### Build Manifest tables####
## define filenames
metadata_types <- c("individual_human", "biospecimen", "assay_rnaSeq", "assay_snpArray","manifest")
meta_files <- paste0("psychENCODE_MDD_", metadata_types,
                     c(rep(".csv", length(metadata_types)-1), ".tsv"))

#Add annotation to metadata files
meta_files <- c(meta_files, annotation_fn)
names(meta_files) <- c(metadata_types, "snpArray_annotation")

## build meta manifest
meta_manifest <- data.frame(
    path = paste0("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/synapse/", meta_files),
    metadataType = names(meta_files),
    fileFormat = ss(meta_files, "\\.", 2)
) %>% mutate(dataSubtype = "metadata",
             BrNum = NA,RNum = NA,
             assay = NA,
             isMultiIndividual = NA,
             isMultiSpecimen = NA)

## build rnaSeq manifest
# build manifest
rnaSeq_manifest <- pd_mdd %>%
    select(RNum, BrNum, read1, read2) %>%
    melt(id.vars = c("RNum", "BrNum")) %>%
    select(RNum, BrNum, path = value) %>%
    mutate(metadataType = NA,
           fileFormat = "fastq",
           assay = "rnaSeq",
           dataSubtype = "raw",
           isMultiIndividual = FALSE,
           isMultiSpecimen = FALSE) %>%
    filter(!is.na(path))

dim(mdd_manifest)
# [1] 1234    7
## build snpArray manifest
snpArray_manifest <- data.frame(
    path = paste0("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/", geno_files),
    fileFormat = gsub("^.*\\.", "", geno_files)
) %>% mutate(metadataType = NA,
             dataSubtype = "processed",
             BrNum = NA,RNum = NA,
             assay = "snpArray",
             isMultiIndividual = TRUE,
             isMultiSpecimen = TRUE)

mdd_all_mani <- meta_manifest %>%
    rbind(rnaSeq_manifest) %>%
    rbind(snpArray_manifest)
message("Dim full manifest:")
dim(mdd_all_mani)
# [1] 1270    9

#### Build Metadata tables####
message('**** Build Metadata tables ****')
## Individual
message("\n** Individual **")
indi_md <- build_metadata("template_individual_human.xlsx", lims, "BrNum", BrNum_mdd)
dim(indi_md)
# [1] 317  33
write.csv(indi_md, file = meta_files[["individual_human"]], row.names = FALSE)

## Biospecimen
message("\n** Biospecimen **")
bio_rna_md <- build_metadata("template_biospecimen.xlsx", pd_mdd, "RNum", RNum_mdd)
dim(bio_rna_md)
# [1] 617  12
bio_geno_md <- build_metadata("template_biospecimen_geno.xlsx", brain_sentrix, "BrNum", BrNum_mdd)
dim(bio_geno_md)
# [1] 316  13
## Combine and save
bio_md <- rbind(bio_rna_md, bio_geno_md)
dim(bio_md)
write.csv(bio_md, file = meta_files[["biospecimen"]], row.names = FALSE)

## Assay
message("\n** Assay **")
## rnaSeq
message("rnaSeq Assay")
rnaSeq_assay_md <- build_metadata("template_assay_rnaSeq.xlsx", pd, "RNum", RNum_mdd)
dim(rnaSeq_assay_md)
# [1] 627  19
write.csv(rnaSeq_assay_md, file = meta_files[["assay_rnaSeq"]], row.names = FALSE)

## snpArray
message("snpArray Assay")
snpArray_assay_md <- build_metadata("template_assay_snpArray.xlsx",fam_df_mdd, "genoSample", genoSample_mdd)
dim(snpArray_assay_md)
# [1] 632   8
write.csv(snpArray_assay_md, file = meta_files[["assay_snpArray"]], row.names = FALSE)

## Manifest
message("\n** Manifest **")
mani_md <- build_metadata("template_manifest.xlsx", mdd_all_mani, "path", mdd_all_mani$path)
map_int(mani_md, ~sum(is.na(.x)))
dim(mani_md)
# [1] 1258   16
write.table(mani_md, file = meta_files[["manifest"]], row.names = FALSE, sep = "\t")

# sgejobs::job_single('build_metadata', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript build_metadata.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# [1] "2020-08-28 13:24:47 EDT"
# user  system elapsed
# 16.960  30.214  46.433
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.0.2 Patched (2020-06-24 r78746)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2020-08-28
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date       lib source
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.0.0)
# Biobase              * 2.48.0   2020-04-27 [2] Bioconductor
# BiocGenerics         * 0.34.0   2020-04-27 [2] Bioconductor
# bitops                 1.0-6    2013-08-17 [2] CRAN (R 4.0.0)
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.0.0)
# cli                    2.0.2    2020-02-28 [1] CRAN (R 4.0.0)
# crayon                 1.3.4    2017-09-16 [2] CRAN (R 4.0.0)
# DelayedArray         * 0.14.0   2020-04-27 [2] Bioconductor
# dplyr                * 0.8.5    2020-03-07 [1] CRAN (R 4.0.0)
# ellipsis               0.3.1    2020-05-15 [2] CRAN (R 4.0.0)
# fansi                  0.4.1    2020-01-08 [2] CRAN (R 4.0.0)
# GenomeInfoDb         * 1.24.0   2020-04-27 [2] Bioconductor
# GenomeInfoDbData       1.2.3    2020-05-18 [2] Bioconductor
# GenomicRanges        * 1.40.0   2020-04-27 [2] Bioconductor
# glue                   1.4.1    2020-05-13 [1] CRAN (R 4.0.0)
# googledrive            1.0.1    2020-05-05 [1] CRAN (R 4.0.0)
# IRanges              * 2.22.1   2020-04-28 [2] Bioconductor
# jaffelab             * 0.99.30  2020-05-21 [1] Github (LieberInstitute/jaffelab@42637ff)
# lattice                0.20-41  2020-04-02 [3] CRAN (R 4.0.2)
# lifecycle              0.2.0    2020-03-06 [1] CRAN (R 4.0.0)
# limma                  3.44.1   2020-04-28 [2] Bioconductor
# magrittr               1.5      2014-11-22 [2] CRAN (R 4.0.0)
# Matrix                 1.2-18   2019-11-27 [3] CRAN (R 4.0.2)
# matrixStats          * 0.56.0   2020-03-13 [2] CRAN (R 4.0.0)
# pillar                 1.4.6    2020-07-10 [1] CRAN (R 4.0.2)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.0.0)
# plyr                   1.8.6    2020-03-03 [2] CRAN (R 4.0.0)
# purrr                * 0.3.4    2020-04-17 [1] CRAN (R 4.0.0)
# R6                     2.4.1    2019-11-12 [2] CRAN (R 4.0.0)
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.0.0)
# RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.0.0)
# Rcpp                   1.0.5    2020-07-06 [1] CRAN (R 4.0.2)
# RCurl                  1.98-1.2 2020-04-18 [2] CRAN (R 4.0.0)
# readxl               * 1.3.1    2019-03-13 [2] CRAN (R 4.0.0)
# reshape2             * 1.4.4    2020-04-09 [2] CRAN (R 4.0.0)
# rlang                  0.4.7    2020-07-09 [1] CRAN (R 4.0.2)
# S4Vectors            * 0.26.1   2020-05-16 [2] Bioconductor
# segmented              1.1-0    2019-12-10 [2] CRAN (R 4.0.0)
# sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.0.0)
# stringi                1.4.6    2020-02-17 [2] CRAN (R 4.0.0)
# stringr                1.4.0    2019-02-10 [2] CRAN (R 4.0.0)
# SummarizedExperiment * 1.18.1   2020-04-30 [2] Bioconductor
# tibble                 3.0.3    2020-07-10 [1] CRAN (R 4.0.2)
# tidyselect             1.1.0    2020-05-11 [2] CRAN (R 4.0.0)
# vctrs                  0.3.2    2020-07-15 [1] CRAN (R 4.0.2)
# withr                  2.2.0    2020-04-20 [1] CRAN (R 4.0.0)
# XVector                0.28.0   2020-04-27 [2] Bioconductor
# zlibbioc               1.34.0   2020-04-27 [2] Bioconductor
#
# [1] /users/lhuuki/R/4.0
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0/R/4.0/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0/R/4.0/lib64/R/library
