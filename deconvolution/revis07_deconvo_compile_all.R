# Louise Huuki-Myers Aug 22
# compile and plot all deconvolution results for MDDseq revisions

library("tidyverse")
library("here")
library("SummarizedExperiment")
library("sessioninfo")

#### Read data ####

load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
pd <- as.data.frame(colData(rse_gene))
pd2 <- pd[,c("RNum", "BrNum", "Experiment", "BrainRegion", "PrimaryDx")]


est_prop <- list()

## Original results - Tran Broad, Bisque
load(here("deconvolution","data","est_prop_Bisque.Rdata"),verbose = TRUE)

est_prop$tran_broad_bisque <- est_prop_bisque$Est.prop.long |>
  separate(Sample, into = c("RNum", "Experiment"), extra = "merge", remove = FALSE) |>
  left_join(pd2) |>
  filter(PrimaryDx != "Bipolar")

## Revis results - Tran Broad, Bisque + hspe



## Revis results - BICCN data, hspe

