library("SummarizedExperiment")
library("tidyverse")
library("here")
library("readxl")

source(here("eqtl", "code", "utils.R"))

## load expression data ##
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
# load(here("exprs_cutoff", "rse_exon.Rdata"), verbose = TRUE)
# load(here("exprs_cutoff", "rse_jxn.Rdata"), verbose = TRUE)
# load(here("exprs_cutoff", "rse_tx.Rdata"), verbose = TRUE)

rd <- as.data.frame(rowData(rse_gene)) %>% select(gencodeID, Symbol)
pd <- as.data.frame(colData(rse_gene))

## build model
mod <- model.matrix(~ PrimaryDx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, data = pd)
colnames(mod)


## TODO move code from summarize nominal
#### get residual expression by region ####
gene_resid_expres_split <- map(splitit(rse_gene$BrainRegion), ~ get_resid_expres(rse_gene[, .x], mod[.x, ]))

gene_resid_long <- map2(gene_resid_expres_split, names(gene_resid_expres_split),
                        ~.x %>%
                          as.data.frame() %>%
                          rownames_to_column("gencodeID") %>%
                          pivot_longer(!gencodeID, names_to = "genoSample", values_to = "expression") %>%
                          mutate(BrainRegion = .y))

gene_resid_long <- do.call("rbind", gene_resid_long)
