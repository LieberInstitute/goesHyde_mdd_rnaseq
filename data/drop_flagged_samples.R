library("sessioninfo")
library("SummarizedExperiment")

## style this script
#styler::style_file("drop_flagged_samples.R", transformers = styler::tidyverse_style(indent_by = 4))

## load rse data
load("rse_gene_GoesZandi_n1093.rda", verbose = TRUE)
load("rse_exon_GoesZandi_n1093.rda", verbose = TRUE)
load("rse_jxn_GoesZandi_n1093.rda", verbose = TRUE)
load("rse_tx_GoesZandi_n1093.rda", verbose = TRUE)

## Find samples to drop
load("../genotype_data/flagged_samples_mddseq.Rdata", verbose = TRUE)

## phenotype data

pd <- colData(rse_gene)
keep_samples <- !pd$RNum %in% checks$SAMPLE_ID
table(keep_samples)
dim(checks)
length(unique(checks$SAMPLE_ID))

#### subset samples
dim(rse_gene)
rse_gene <- rse_gene[, keep_samples]
dim(rse_gene)

rse_exon <- rse_exon[, keep_samples]
rse_jxn <- rse_jxn[, keep_samples]
rse_tx <- rse_tx[, keep_samples]

## save data
dir.create("without_flagged", showWarnings = FALSE)
save(rse_gene, file = file.path("without_flagged", "rse_gene_bp_mdd.Rdata"))
save(rse_exon, file = file.path("without_flagged", "rse_exon_bp_mdd.Rdata"))
save(rse_jxn, file = file.path("without_flagged", "rse_jxn_bp_mdd.Rdata"))
save(rse_tx, file = file.path("without_flagged", "rse_tx_bp_mdd.Rdata"))


# sgejobs::job_single('drop_flagged_samples', create_shell = TRUE, memory = '40G', command = "Rscript drop_flagged_samples.R")


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
