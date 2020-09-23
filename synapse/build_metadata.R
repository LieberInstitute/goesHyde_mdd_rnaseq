library("readxl")
library("sessioninfo")
library("SummarizedExperiment")
library("jaffelab")
library("dplyr")
library("reshape2")
library("purrr")
libary("here")

## style this script
# styler::style_file("build_metadata.R", transformers = styler::tidyverse_style(indent_by = 4))

build_metadata <- function(template_xlsx, data, id_col, id) {
    dict <- read_excel(template_xlsx, sheet = "Dictionary")
    dict_id <- dict$key[match(id_col, dict$col)]
    message(paste("data ID:", id_col, "template id:", dict_id))
    # get extract variable cols from data
    data_hasCol <- !is.na(dict$col) & dict$col != "?"
    data_col <- dict$col[data_hasCol]
    
    # message(paste(colnames(data), collapse = ","))
    # message(paste(data_col, collapse = ","))
    # check colnames match
    col_match <- data_col %in% colnames(data)
    if (!all(col_match)) {
        missing <- data_col[!col_match]
        message("Missing cols:", paste(missing, collapse = ","))
        return(NULL)
    } else {
        message(sum(col_match), " matches in data")
    }
    
    # build data Variable
    dataV <- data[data[[id_col]] %in% id, ]
    dataV <- dataV[, data_col]
    colnames(dataV) <- dict$key[data_hasCol]
    dataV <- replace_values(template_xlsx, dataV)
    # build data Same
    dataS <- t(data.frame(dict$value[!data_hasCol]))
    dataS <- do.call("rbind", replicate(length(id), dataS, simplify = FALSE))
    dataS <- cbind(dataS, id)
    colnames(dataS) <- c(dict$key[!data_hasCol], dict_id)
    # build data All
    dataA <- merge(dataV, dataS, by = dict_id)
    temp <- read_excel(template_xlsx, sheet = "Template")
    
    meta_data <- rbind(temp, dataA)
    return(meta_data)
}

get_fastq_info <- function(fastq) {
    l <- system(paste0("zcat ", fastq, ' | grep "@" | head -n 1'), intern = TRUE) %>%
        strsplit(" ") %>%
        unlist()
    
    l1 <- strsplit(l, ":") %>% unlist()
    info_names <- c(
        "instrument", "rna_id", "flow_cell", "flowcell_lane",
        "title_number", "x_cord", "y_cord", "pair", "filtered",
        "control_bits", "index_seq"
    )
    names(l1) <- info_names
    return(l1)
}

replace_value <- function(value_row, dataV) {
    cn <- value_row[1]
    v <- value_row[2]
    lv <- value_row[3]
    dataV[[cn]][dataV[[cn]] == lv] <- v
    return(dataV)
}

make_value_df <- function(template_xlsx) {
    value_df <- read_excel(template_xlsx, "Values") %>%
        select(key, value, LIBD_value) %>%
        filter(
            !is.na(LIBD_value),
            value != LIBD_value
        ) %>%
        as.data.frame()
    message(nrow(value_df), " values to replace")
    return(value_df)
}

replace_values <- function(template_xlsx, dataV) {
    value_df <- make_value_df(template_xlsx)
    nr <- nrow(value_df)
    if (nr == 0) {
        return(dataV)
    } else {
        for (row in 1:nrow(value_df)) {
            dataV <- replace_value(unlist(value_df[row, ]), dataV)
        }
        return(dataV)
    }
}

#### load data ####
lims <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/shiny070220.csv")
pd <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/pd_swap.csv", as.is = TRUE)

# add brain weight to pd
brain_weight <- lims %>% select(BrNum, `Brain.Weight..gram.`)

pd_mdd <- read.csv(here("data","raw_GoesZandi_pd.csv")) %>%
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

# build manifest
mdd_manifest <- pd_mdd %>%
    select(RNum, BrNum, read1, read2) %>%
    melt(id.vars = c("RNum", "BrNum")) %>%
    select(RNum, BrNum, path = value) %>%
    mutate(metadataType = NA, 
           fileFormat = "fastq",
           assay = "rnaSeq",
           dataSubtype = "raw") %>%
    filter(!is.na(path))

# fastq info
fastq_info <- map(mdd_manifest$path, get_fastq_info)
fastq_info_df <- fastq_info %>%
    as.data.frame() %>%
    t()
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

#### Build metadata ####
#define filenames
indi_md_fn <- "psychENCODE_MDD_individual_human.csv"
bio_md_fn <- "psychENCODE_MDD_biospecimen.csv"
assay_md_fn <- "psychENCODE_MDD_assay_rnaSeq.csv"
mani_md_fn <- "psychENCODE_MDD_manifest.tsv"

# Individual
message("\nIndividual")
indi_md <- build_metadata("template_individual_human.xlsx", lims, "BrNum", BrNum_mdd)
dim(indi_md)
# [1] 317  33
write.csv(indi_md, file = indi_md_fn, row.names = FALSE)

# Biospecimen
message("\nBiospecimen")
bio_md <- build_metadata("template_biospecimen.xlsx", pd_mdd, "RNum", RNum_mdd)
dim(bio_md)
# [1] 627  12
write.csv(bio_md, file = bio_md_fn, row.names = FALSE)

# Assay
message("\nAssay")
assay_md <- build_metadata("template_assay_rnaSeq.xlsx", pd, "RNum", RNum_mdd)
dim(assay_md)
# [1] 627  19
write.csv(assay_md, file = assay_md_fn, row.names = FALSE)

# Manifest
message("\nManifest")

meta_files <- c(indi_md_fn, bio_md_fn, assay_md_fn, mani_md_fn)
meta_paths <- paste0("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/synapse/", meta_files)
meta_manifest <- data.frame(
    meta_paths,
    c("individual", "biospecimen", "assay", "manifest"),
    ss(meta_files, "\\.", 2),
    rep("metadata",4)
)
colnames(meta_manifest) <- c("path", "metadataType", "fileFormat","dataSubtype")

# build genotype file manifest
geno_files <- scan("../genotype_data/manifest.txt", what="character", sep="\n")
geno_n <- length(geno_files)
geno_paths <- paste0("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/", geno_files)
geno_manifest <- data.frame(
    geno_paths,
    gsub("^.*\\.", "", geno_files),
    rep("Genotyping",geno_n),
    rep("processed", geno_n)
)
colnames(geno_manifest) <- c("path", "fileFormat", "assay", "dataSubtype")


mdd_all_mani <- meta_manifest %>%
    mutate(assay = NA) %>%
    rbind(geno_manifest %>% mutate(metadataType = NA)) %>%
    mutate(BrNum = NA,RNum = NA) %>%
    rbind(mdd_manifest)

mani_md <- build_metadata("template_manifest.xlsx", mdd_all_mani, "path", mdd_all_mani$path)
map_int(mani_md, ~sum(is.na(.x)))
dim(mani_md)
# [1] 1258   16
write.table(mani_md, file = mani_md_fn, row.names = FALSE, sep = "\t")

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
