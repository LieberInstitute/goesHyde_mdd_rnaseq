library(SummarizedExperiment)
library(eQTLUtils)
library(sessioninfo)
library(Tidyverse)
library(here)

source(here("eqtl","risk_snps","rse_to_bed.R"))

## load
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
# load(here("exprs_cutoff", "rse_exon.Rdata"), verbose = TRUE)
# load(here("exprs_cutoff", "rse_jxn.Rdata"), verbose = TRUE)
# load(here("exprs_cutoff", "rse_tx.Rdata"), verbose = TRUE)

## Split gene data
regions <- c(amyg = "Amygdala", sacc = "sACC")
rse_gene_split <- map(regions, ~rse_gene[,rse_gene$BrainRegion == .x])

#### Covariate Data ####
load(here("eqtl", "genomewide", "rdas", "pcs_4features_amyg.rda"), verbose = TRUE)

map(list(genePCs, exonPCs, jxnPCs, txPCs), function(pc){
  pc <- t(pc)
  table(ss(colnames(pc), "_") == rse_gene_split$amyg$RNum)
})

ss(rownames(genePCs),"_")[(ss(rownames(genePCs),"_") %in% rse_gene$RNum)]


#### Create Phenotype Data ####
fn_gene <- map(regions, ~paste0("pheno_data/gene_",.x,".bed"))

walk2(rse_gene_split, fn_gene, function(rse,fn){
  bed_gene <- rse_to_bed(rse)
  write.table(bed_gene, fn, sep = "\t", 
              quote = FALSE, row.names = FALSE)
})

## Create shell script to zip data
commands <- map(fn_gene, ~paste0("bgzip ", .x," && tabix -p bed ", .x,".gz"))
sgejobs::job_single("bed_bgzip", create_shell = TRUE, queue= 'bluejay', memory = '100G', 
                    command = paste(commands, collapse = '\n'))


