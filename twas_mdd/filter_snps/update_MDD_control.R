
library(sva)
library(lmerTest)
library(SummarizedExperiment)
library(jaffelab)
library(WGCNA)
library(broom)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(readxl)
library(RColorBrewer)
library(sessioninfo)
library(here)

rm(list = ls())


capabilities()

## load data

load("../exprs_cutoff/rse_gene.Rdata", verbose = TRUE)

rse_gene <- rse_gene[, rse_gene$PrimaryDx %in% c("Control", "MDD")]

rse_gene_Amygdala <- rse_gene[, rse_gene$BrainRegion %in% c("Amygdala")]
rse_gene_sACC <- rse_gene[, rse_gene$BrainRegion %in% c("sACC")]

setwd("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/twas_mdd/filter_snps")

write.table(rse_gene_sACC$genoSample, file = "samples_sACC_MDD_control.txt", quote = F, row.names = FALSE, col.names = FALSE)

write.table(rse_gene_Amygdala$genoSample, file = "samples_Amyg_MDD_control.txt", quote = F, row.names = FALSE, col.names = FALSE)
