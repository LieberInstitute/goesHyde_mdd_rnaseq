library(SummarizedExperiment)
library(sessioninfo)
library(tidyverse)
library(jaffelab)
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
# load(here("eqtl", "genomewide", "rdas", "pcs_4features_amyg.rda"), verbose = TRUE)
load(here("eqtl", "pcs_4features_gene_temp.rda"), verbose = TRUE)

corner(genePCs)

all(rownames(genePCs) == colnames(rse_gene))

# map(list(genePCs, exonPCs, jxnPCs, txPCs), function(pc){
covars <- map2(rse_gene_split, regions,  ~map(list(genePCs), function(pc){
  message(.y)
  pc <- pc[colnames(.x),]
  rownames(pc) <- .x$genoSample
  pc <- t(pc)
  pc <- as.data.frame(pc) %>% rownames_to_column("id")
  write.table(pc,
              file= here("eqtl","risk_snps","pheno_data", paste0("covariates_gene_",.y,".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  return(pc)
}))

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

sgejobs::job_single("tensorQTL_test", create_shell = TRUE, queue= 'bluejay', memory = '10G', 
                    command = 
                      "python3 -m tensorqtl LIBD_maf01_gwas_BPD pheno_data/gene_Amygdala.bed.gz gene_Amygdala \
    --covariates pheno_data/covariates_gene_Amygdala.txt \
    --mode cis")



