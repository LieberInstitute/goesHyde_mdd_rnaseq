
library(jaffelab)
library(SummarizedExperiment)
library(tidyverse)
library(reshape2)
library(GGally)
library(sessioninfo)
library(here)

source(here("main_colors.R"))

## load rse
load(here("genotype_data","goesHyde_bipolarMdd_Genotypes_mds.rda"), verbose = TRUE)

load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)
pd <- as.data.frame(colData(rse_gene))
snpPC <- pd[,grepl( "snpPC", colnames(pd))]

## Just PCs
pair_plot <- ggpairs(mds, aes(alpha = 0.4))
ggsave(pair_plot, filename = "plots/snpPC_explore.png")

## with Dx
snpPC <- cbind(pd[,"PrimaryDx",FALSE], snpPC)
pair_plot <- ggpairs(snpPC, aes(color = PrimaryDx, alpha = 0.4)) 
# +
#   scale_color_manual(values = mdd_Dx_colors) +
#   scale_fill_manual(values = mdd_Dx_colors)

ggsave(pair_plot, filename = "plots/snpPC_explore_Dx.png", width = 10, height = 10)

#sgejobs::job_single('snpPC_explore', create_shell = TRUE, memory = '5G', command = "Rscript snpPC_explore.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

