### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)
library(sessioninfo)
library(here)

## load
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_exon.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_jxn.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_tx.Rdata"), verbose = TRUE)

pd <- colData(rse_gene)

load(here("genotype_data","goesHyde_bipolarMdd_Genotypes_mds.rda"), verbose = TRUE)

#### Model ####
message(Sys.time(), " Get statistical model")
pd$PrimaryDx <- factor(pd$PrimaryDx,
                       levels = c("Control", "Bipolar", "MDD")
)

mod <- model.matrix(~ PrimaryDx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, data = pd)
colnames(mod)

#### calculate rpkm ####
message(Sys.time(), " Get rpkm values")
geneRpkm <- recount::getRPKM(rse_gene, "Length")
rm(rse_gene)

exonRpkm <- recount::getRPKM(rse_exon, "Length")
rm(rse_exon)

rowData(rse_jxn)$Length <- 100
jxnRp10m <- recount::getRPKM(rse_jxn, "Length")
rm(rse_jxn)

txTpm <- assays(rse_tx)$tpm
rm(rse_tx

#### do PCA ####
message(Sys.time(), " Do PCA")

pca_rda_file <- here("eqtl", "pcs_4features.rda")

pcaGene <- prcomp(t(log2(geneRpkm + 1)))
kGene <- num.sv(log2(geneRpkm + 1), mod)
genePCs <- pcaGene$x[, 1:kGene]

pcaExon <- prcomp(t(log2(exonRpkm + 1)))
kExon <- num.sv(log2(exonRpkm + 1), mod, vfilter = 50000)
exonPCs <- pcaExon$x[, 1:kExon]

pcaJxn <- prcomp(t(log2(jxnRp10m + 1)))
kJxn <- num.sv(log2(jxnRp10m + 1), mod, vfilter = 50000)
jxnPCs <- pcaJxn$x[, 1:kJxn]

pcaTx <- prcomp(t(log2(txTpm + 1)))
kTx <- num.sv(log2(txTpm + 1), mod, vfilter = 50000)
txPCs <- pcaTx$x[, 1:kTx]

save(genePCs, exonPCs, jxnPCs, txPCs, file = pca_rda_file)
  
# sgejobs::job_single("get_pcs_4features", memory = "150G",create_shell = TRUE, command = "Rscript get_pcs_4features.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

