### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)
library(purrr)
library(sessioninfo)
library(here)

## load
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_exon.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_jxn.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_tx.Rdata"), verbose = TRUE)

pd <- colData(rse_gene)

# load(here("genotype_data","goesHyde_bipolarMdd_Genotypes_mds.rda"), verbose = TRUE)

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

exonRpkm <- recount::getRPKM(rse_exon, "Length")
rm(rse_exon)

rowData(rse_jxn)$Length <- 100
jxnRp10m <- recount::getRPKM(rse_jxn, "Length")
rm(rse_jxn)

txTpm <- assays(rse_tx)$tpm
rm(rse_tx)

## Set up map
region_samples <- splitit(rse_gene$BrainRegion)
map(region_samples, length)

map(region_samples, ~ dim(rse_gene[, .x]))

#### do PCA ####
(pca_rda_dir <- here("eqtl", "data", "featuresPCs/"))

## genePCs
# genePCs <- map2(region_samples, names(region_samples), function(r_samples, r_name){
#   message(Sys.time(), paste(" gene", r_name, "PCA"))
#
#   rpkm_region <- geneRpkm[,r_samples]
#   mod_region <- mod[r_samples,]
#   pca <- prcomp(t(log2(rpkm_region + 1)))
#   message("start num sv")
#   k <- num.sv(log2(rpkm_region + 1), mod_region)
#   PCs <- pca$x[, 1:k]
#
#   save(PCs, file = paste0(pca_rda_dir, "genePCs",r_name,".rda"))
#   return(PCs)
# })

## Other Featuers w/ vfilter
otherPCs <- map2(list(exon = exonRpkm, jxn = jxnRp10m, tx = txTpm), c("exon", "jxn", "tx"), function(f_rpkm, f_name) {
    map2(region_samples, names(region_samples), function(r_samples, r_name) {
        message(Sys.time(), paste("", f_name, r_name, "PCA"))

        rpkm_region <- f_rpkm[, r_samples]
        mod_region <- mod[r_samples, ]
        pca <- prcomp(t(log2(rpkm_region + 1)))
        k <- num.sv(log2(rpkm_region + 1), mod_region, vfilter = 50000)
        PCs <- pca$x[, 1:k]

        save(PCs, file = paste0(pca_rda_dir, f_name, "PCs", r_name, ".rda"))
        return(PCs)
    })
})

save(otherPCs, file = here("eqtl", "data", "featuresPCs", "pcs_4features.Rdata"))
# save(genePCs, otherPCs, file = here("eqtl", "data", "featuresPCs","pcs_4features.Rdata"))
# sgejobs::job_single("get_pcs_4features", memory = "150G",create_shell = TRUE, command = "Rscript get_pcs_4features.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
