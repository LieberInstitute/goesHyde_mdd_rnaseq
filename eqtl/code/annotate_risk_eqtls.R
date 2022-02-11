library(SummarizedExperiment)
library(tidyverse)
library(here)

load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_exon.Rdata"), verbose = TRUE)
# load(here("exprs_cutoff", "rse_jxn.Rdata"), verbose = TRUE)
# load(here("exprs_cutoff", "rse_tx.Rdata"), verbose = TRUE)

## Risk eQTL files
fn <- list.files(path = here("eqtl", "data", "risk_snps_eqtl"), full.names = TRUE)
risk_eqtl <- read.csv(fn[[1]])

head(risk_eqtl)

rowData(rse_gene)

rd <- rowData(rse_exon) %>%
  as.data.frame() %>%
  select(gencodeID, Symbol, gene_type, Class) %>%
  rownames_to_column("phenotype_id")
  
head(rd)

risk_eqtl_anno <- left_join(risk_eqtl, rd)
head(risk_eqtl_anno)

risk_eqtl_anno %>% 
  filter(FDR < 0.05) %>% 
  count(Symbol) %>%
  arrange(-n)
