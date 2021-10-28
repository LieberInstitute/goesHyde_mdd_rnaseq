### libraries
library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(segmented)
library(here)

## load
load(here("exprs_cutoff","rse_gene.Rdata"), verbose = TRUE)
load(here("exprs_cutoff","rse_exon.Rdata"), verbose = TRUE)
load(here("exprs_cutoff","rse_jxn.Rdata"), verbose = TRUE)
load(here("exprs_cutoff","rse_tx.Rdata"), verbose = TRUE)

load(here("genotype_data","goesHyde_bipolarMdd_Genotypes_mds.rda"), verbose = TRUE)
## Check that all BrNum are included
message("# BrNum in rse: ", length(unique(rse_gene$BrNum)),
        ", # BrNum in mds: ", nrow(mds))
message("Uncommon BrNum: ", length(setdiff(rownames(mds), unique(rse_gene$BrNum))))

## reorders according to rse_gene$BrNum  (matched the order of rse_gene)
mds = mds[rse_gene$BrNum,1:5]
dim(mds)

all(rse_gene$BrNum == ss(rownames(mds),"\\."))

colData(rse_gene)[,paste0("snpPC", 1:5)] <- mds
colData(rse_exon)[,paste0("snpPC", 1:5)] <- mds
colData(rse_jxn)[,paste0("snpPC", 1:5)] <- mds
colData(rse_tx)[,paste0("snpPC", 1:5)] <- mds

save(rse_gene, file = here('exprs_cutoff','rse_gene.Rdata'))
save(rse_exon, file = here('exprs_cutoff','rse_exon.Rdata'))
save(rse_jxn, file = here('exprs_cutoff','rse_jxn.Rdata'))
save(rse_tx, file = here('exprs_cutoff','rse_tx.Rdata'))
