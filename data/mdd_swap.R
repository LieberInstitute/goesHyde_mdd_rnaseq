library(jaffelab)
library(VariantAnnotation)
library(readxl)
library(Rsamtools)
library(dplyr)
library(SummarizedExperiment)
library(ggplot2)
library(sessioninfo)
#### Load Data ####
load("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/corLong2.Rdata", verbose = TRUE)
pd <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/pd_swap.csv", as.is=TRUE)

# MDD data
load("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/data/rse_gene_raw_GoesZandi_n1174.rda", verbose = TRUE)
mdd_pd <- colData(rse_gene)
dim(mdd_pd) #[1] 1174   13
length(unique(mdd_pd$BrNum)) #[1] 607

table(rse_gene$Experiment)
# GoesMDD ZandiBPD 
# 634      540

#### Swap and drop rna samples ###
# load pd data
pd <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/pd_swap.csv") %>%
  rename(AgeDeath = Age, BrainRegion = Region, PrimaryDx = Dx, Experiment = Dataset) %>%
  select(colnames(mdd_pd)) %>%
  filter(RNum %in% mdd_pd$RNum) 

pd$Experiment <- factor(pd$Experiment, levels = c("psychENCODE_MDD","psychENCODE_BP"))
table(pd$Experiment)
# psychENCODE_MDD  psychENCODE_BP 
# 634             540 

# drop 10 overlapping BPD Control samples
pd <- pd %>%
  group_by(RNum) %>% 
  arrange(Experiment) %>%
  slice(1) %>% 
  ungroup()

pd %>% count(Experiment, BrNum == "drop")
# Experiment      `BrNum == "drop"`     n
# <fct>           <lgl>             <int>
#   1 psychENCODE_MDD FALSE               631
# 2 psychENCODE_MDD TRUE                  3
# 3 psychENCODE_BP  FALSE               523
# 4 psychENCODE_BP  TRUE                  7


#### Match with Brain_Sentrix samples ####
#brain sentrix info
brain_sentrix<- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/brain_sentrix_swap.csv") %>%
  filter(BrNum != "drop", BrNum %in% pd$BrNum)

#load batch priority (lower is better)
array_priority <- read.csv("/dcl01/ajaffe/data/lab/brain_swap/DNA_genotyping_array_priorities.csv")
array_priority <- array_priority[,c(1,4)]

# swap filter prioritize
brain_sentrix1 <- brain_sentrix %>%
  left_join(array_priority, by = "Batch") %>%
  group_by(BrNum) %>%
  mutate(nBr = n(),
         n_batch = length(unique(Batch)),
         MinPriority = min(Priority)) %>%
  filter(Priority == MinPriority) %>%
  arrange(ID)%>%
  slice(1) %>%
  ungroup() %>%
  select(genoSample = ID, BrNum, Batch)

# should only have 1 row per BrNum
nrow(brain_sentrix1) ==length(unique(brain_sentrix1$BrNum)) 

#### Add brain sentrix ID to mdd_data ####
pd <- pd %>% left_join(brain_sentrix1, by = "BrNum")

# check for good cor in corLong2
corLong2_mdd <- corLong2 %>% select(RNum, genoSample, cor) %>%
  filter(genoSample %in% pd$genoSample,
         RNum %in% pd$RNum) %>%
  group_by(RNum,genoSample) %>%
  arrange(-cor)%>%
  slice(1) %>%
  ungroup

pd_check <- pd %>% left_join(corLong2_mdd, by = c("RNum","genoSample"))

table(pd_check$cor >= 0.59 & !is.na(pd_check$cor))
#Drop 18 samples 
# FALSE  TRUE 
# 18  1146 

# drop break down by Experiment
pd_check %>%
  count(Experiment, dna_match = cor >= 0.59 & !is.na(cor), good_rna = BrNum != "drop")
# Experiment dna_match good_rna     n
# <chr>      <lgl>     <lgl>    <int>
# 1 GoesMDD    FALSE     FALSE        3
# 2 GoesMDD    FALSE     TRUE         4
# 3 GoesMDD    TRUE      TRUE       627
# 4 ZandiBPD   FALSE     FALSE        7
# 5 ZandiBPD   FALSE     TRUE         4
# 6 ZandiBPD   TRUE      TRUE       519

pd_good <- pd_check %>% 
  filter(BrNum != "drop" & cor >=0.59) 

#filter out brains not in lims
lims <- read.csv("/dcl01/lieber/RNAseq/Datasets/BrainGenotyping_2018/SampleFiles/preBrainStorm/shiny070220.csv")
br_missing_lims <- pd_good$BrNum[!pd_good$BrNum %in% lims$BrNum] %>% unique
message(paste("Brains missing from lims:",length(br_missing_lims)))
message(paste(br_missing_lims, collapse = " ,"))

pd_good <- pd_good %>% filter(BrNum %in% lims$BrNum)

# create final pd table
pd <- pd_good %>% select(-cor)

dna_rna_cor_histo <- pd_good %>%
  ggplot(aes(cor)) +
  geom_histogram() +
  labs(title = "DNA vs. RNA cor",
       subtitle = paste0(nrow(pd)," good samples for MDD"))

ggsave(filename = "MDDsamples_cor_histo.jpg", plot = dna_rna_cor_histo)



#MDD fastq files
fastq_mdd <- paste0("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/FASTQ_merged/",pd$RNum[pd$Experiment == "psychENCODE_MDD"],".fastq.gz")
fastq_mdd2 <- paste0("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/FASTQ_merged/",pd$RNum[pd$Experiment == "psychENCODE_MDD"],"_read2.fastq.gz")
#check valid paths
fe_mdd <- file.exists(c(fastq_mdd, fastq_mdd2))
table(fe_mdd)
# TRUE 
# 1254 

fastd_df <- data.frame(pd$RNum[pd$Experiment == "psychENCODE_MDD"], fastq_mdd, fastq_mdd2)
colnames(fastd_df) <- c("RNum", "read1", "read2")

# Add fastq paths 
pd <- pd %>% 
  left_join(fastd_df) %>%
  arrange(BrNum) %>%
  arrange(Experiment)

n = nrow(pd) 

fn = paste0("GoesMDD_pd_n",n,".csv")
write.csv(pd, fn)

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# sgejobs::job_single('mdd_swap', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript mdd_swap.R")

# [1] "Reproducibility information:"
# [1] "2020-08-25 16:54:40 EDT"
# user  system elapsed 
# 79.346   2.683  83.447 
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
# date     2020-08-25                                 
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date       lib source                                   
# AnnotationDbi          1.50.0   2020-04-27 [2] Bioconductor                             
# askpass                1.1      2019-01-13 [2] CRAN (R 4.0.0)                           
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.0.0)                           
# Biobase              * 2.48.0   2020-04-27 [2] Bioconductor                             
# BiocFileCache          1.12.0   2020-04-27 [2] Bioconductor                             
# BiocGenerics         * 0.34.0   2020-04-27 [2] Bioconductor                             
# BiocParallel           1.22.0   2020-04-27 [2] Bioconductor                             
# biomaRt                2.44.0   2020-04-27 [2] Bioconductor                             
# Biostrings           * 2.56.0   2020-04-27 [2] Bioconductor                             
# bit                    1.1-15.2 2020-02-10 [2] CRAN (R 4.0.0)                           
# bit64                  0.9-7    2017-05-08 [2] CRAN (R 4.0.0)                           
# bitops                 1.0-6    2013-08-17 [2] CRAN (R 4.0.0)                           
# blob                   1.2.1    2020-01-20 [2] CRAN (R 4.0.0)                           
# BSgenome               1.56.0   2020-04-27 [2] Bioconductor                             
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.0.0)                           
# cli                    2.0.2    2020-02-28 [1] CRAN (R 4.0.0)                           
# colorspace             1.4-1    2019-03-18 [2] CRAN (R 4.0.0)                           
# crayon                 1.3.4    2017-09-16 [2] CRAN (R 4.0.0)                           
# curl                   4.3      2019-12-02 [2] CRAN (R 4.0.0)                           
# DBI                    1.1.0    2019-12-15 [2] CRAN (R 4.0.0)                           
# dbplyr                 1.4.3    2020-04-19 [2] CRAN (R 4.0.0)                           
# DelayedArray         * 0.14.0   2020-04-27 [2] Bioconductor                             
# digest                 0.6.25   2020-02-23 [1] CRAN (R 4.0.0)                           
# dplyr                * 0.8.5    2020-03-07 [1] CRAN (R 4.0.0)                           
# ellipsis               0.3.1    2020-05-15 [2] CRAN (R 4.0.0)                           
# fansi                  0.4.1    2020-01-08 [2] CRAN (R 4.0.0)                           
# farver                 2.0.3    2020-01-16 [2] CRAN (R 4.0.0)                           
# GenomeInfoDb         * 1.24.0   2020-04-27 [2] Bioconductor                             
# GenomeInfoDbData       1.2.3    2020-05-18 [2] Bioconductor                             
# GenomicAlignments      1.24.0   2020-04-27 [2] Bioconductor                             
# GenomicFeatures        1.40.0   2020-04-27 [2] Bioconductor                             
# GenomicRanges        * 1.40.0   2020-04-27 [2] Bioconductor                             
# ggplot2              * 3.3.0    2020-03-05 [2] CRAN (R 4.0.0)                           
# glue                   1.4.1    2020-05-13 [1] CRAN (R 4.0.0)                           
# googledrive            1.0.1    2020-05-05 [1] CRAN (R 4.0.0)                           
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.0.0)                           
# hms                    0.5.3    2020-01-08 [2] CRAN (R 4.0.0)                           
# httr                   1.4.2    2020-07-20 [1] CRAN (R 4.0.2)                           
# IRanges              * 2.22.1   2020-04-28 [2] Bioconductor                             
# jaffelab             * 0.99.30  2020-05-21 [1] Github (LieberInstitute/jaffelab@42637ff)
# labeling               0.3      2014-08-23 [2] CRAN (R 4.0.0)                           
# lattice                0.20-41  2020-04-02 [3] CRAN (R 4.0.2)                           
# lifecycle              0.2.0    2020-03-06 [1] CRAN (R 4.0.0)                           
# limma                  3.44.1   2020-04-28 [2] Bioconductor                             
# magrittr               1.5      2014-11-22 [2] CRAN (R 4.0.0)                           
# Matrix                 1.2-18   2019-11-27 [3] CRAN (R 4.0.2)                           
# matrixStats          * 0.56.0   2020-03-13 [2] CRAN (R 4.0.0)                           
# memoise                1.1.0    2017-04-21 [2] CRAN (R 4.0.0)                           
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.0.0)                           
# openssl                1.4.2    2020-06-27 [1] CRAN (R 4.0.2)                           
# pillar                 1.4.6    2020-07-10 [1] CRAN (R 4.0.2)                           
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.0.0)                           
# prettyunits            1.1.1    2020-01-24 [2] CRAN (R 4.0.0)                           
# progress               1.2.2    2019-05-16 [2] CRAN (R 4.0.0)                           
# purrr                  0.3.4    2020-04-17 [1] CRAN (R 4.0.0)                           
# R6                     2.4.1    2019-11-12 [2] CRAN (R 4.0.0)                           
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.0.0)                           
# rappdirs               0.3.1    2016-03-28 [2] CRAN (R 4.0.0)                           
# RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.0.0)                           
# Rcpp                   1.0.5    2020-07-06 [1] CRAN (R 4.0.2)                           
# RCurl                  1.98-1.2 2020-04-18 [2] CRAN (R 4.0.0)                           
# readxl               * 1.3.1    2019-03-13 [2] CRAN (R 4.0.0)                           
# rlang                  0.4.7    2020-07-09 [1] CRAN (R 4.0.2)                           
# Rsamtools            * 2.4.0    2020-04-27 [2] Bioconductor                             
# RSQLite                2.2.0    2020-01-07 [2] CRAN (R 4.0.0)                           
# rtracklayer            1.48.0   2020-04-27 [2] Bioconductor                             
# S4Vectors            * 0.26.1   2020-05-16 [2] Bioconductor                             
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.0.0)                           
# segmented              1.1-0    2019-12-10 [2] CRAN (R 4.0.0)                           
# sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.0.0)                           
# stringi                1.4.6    2020-02-17 [2] CRAN (R 4.0.0)                           
# stringr                1.4.0    2019-02-10 [2] CRAN (R 4.0.0)                           
# SummarizedExperiment * 1.18.1   2020-04-30 [2] Bioconductor                             
# tibble                 3.0.3    2020-07-10 [1] CRAN (R 4.0.2)                           
# tidyselect             1.1.0    2020-05-11 [2] CRAN (R 4.0.0)                           
# utf8                   1.1.4    2018-05-24 [2] CRAN (R 4.0.0)                           
# VariantAnnotation    * 1.34.0   2020-04-27 [2] Bioconductor                             
# vctrs                  0.3.2    2020-07-15 [1] CRAN (R 4.0.2)                           
# withr                  2.2.0    2020-04-20 [1] CRAN (R 4.0.0)                           
# XML                    3.99-0.3 2020-01-20 [2] CRAN (R 4.0.0)                           
# XVector              * 0.28.0   2020-04-27 [2] Bioconductor                             
# zlibbioc               1.34.0   2020-04-27 [2] Bioconductor                             
# 
# [1] /users/lhuuki/R/4.0
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0/R/4.0/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0/R/4.0/lib64/R/library