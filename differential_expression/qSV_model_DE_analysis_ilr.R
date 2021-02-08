#####################################
# from /dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/R_scripts/DEG_scripts/qSV_model_analysis_DLPFC_PTSD.R
#####################################

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)
library(edgeR)
library(purrr)
library(sessioninfo)
library(here)

source("run_DE.R")

##### Load Data ####

## load rse
load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)

pd = colData(rse_gene)

## qSV data
load(here("differential_expression","data","qSV_mat.Rdata"), verbose = TRUE)

## load ilr data
load(here("deconvolution","data","est_prop_ilr.Rdata"), verbose = TRUE)

#### Define models ####
regions <- list(sacc = "sACC", amyg = "Amygdala")

modSep <- map(regions, ~cbind(model.matrix(~PrimaryDx + AgeDeath + Sex  +
                                      snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
                                      mitoRate + rRNA_rate +totalAssignedGene + RIN + ERCCsumLogErr,
                                    data= pd[pd$BrainRegion == .x,]),
                              qSV_mat[pd$BrainRegion == .x,]))

## Add ilr terms
modSep_ilr <- map2(modSep, est_prop_ilr[c("sacc_specific","amyg_specific")], ~cbind(.x, .y))

map(modSep_ilr, colnames)

## save
save(modSep, modSep_ilr, file = "differential_models_ilr.Rdata")

#### Gene ####
message("\nGENE")
rse_gene_sep <- map(regions, ~rse_gene[,rse_gene$BrainRegion == .x])
map(rse_gene_sep, dim)

outGene <- map2(rse_gene_sep, modSep, run_DE)
outGene_ilr <- map2(rse_gene_sep, modSep_ilr, run_DE)

save(outGene, outGene_ilr, file = here("differential_expression","data", "qSVA_MDD_gene_DEresults_ilr.rda"))


#### Exon ####
message("\nEXON")
load(here('exprs_cutoff','rse_exon.Rdata'), verbose=TRUE)
rse_exon <- map(regions, ~rse_exon[,rse_exon$BrainRegion == .x])
map(rse_exon, dim)

outExon <- map2(rse_exon, modSep, run_DE)
save(outExon, file = here("differential_expression","data", "qSVA_MDD_exon_DEresults.rda"))
rm(rse_exon)

#### Jxn ####
message("\nJXN")
load(here('exprs_cutoff','rse_jxn.Rdata'), verbose=TRUE)
rse_jxn <- map(regions, ~rse_jxn[,rse_jxn$BrainRegion == .x])
map(rse_jxn, dim)

outjxn <- map2(rse_jxn, modSep, run_DE)
save(outjxn, file = here("differential_expression","data", "qSVA_MDD_jxn_DEresults.rda"))
rm(rse_jxn)

#### Tx ####
message("\nTX")
load(here('exprs_cutoff','rse_tx.Rdata'), verbose=TRUE)
rse_tx <- map(regions, ~rse_tx[,rse_tx$BrainRegion == .x])
map(rse_tx, dim)

outtx <- map2(rse_tx, modSep, ~run_DE(.x, .y, run_voom = FALSE))
save(outtx, file = here("differential_expression","data", "qSVA_MDD_tx_DEresults.rda"))


#sgejobs::job_single('qSV_model_DE_analysis_ilr', create_shell = TRUE, memory = '80G', command = "Rscript qSV_model_DE_analysis_ilr.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
