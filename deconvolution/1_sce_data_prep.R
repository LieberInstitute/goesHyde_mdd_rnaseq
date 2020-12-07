library(SingleCellExperiment)
library(jaffelab)
library(here)
library(GenomicFeatures)

#### Load and filter data ####
## Load rse_gene data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sACC data
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose = TRUE)
sce.sacc$uniqueID <- paste0(sce.sacc$donor, "_", sce.sacc$Barcode)
colnames(sce.sacc) <- sce.sacc$uniqueID

sce.sacc <- sce.sacc[,sce.sacc$cellType != "Ambig.lowNtrxts",]
sce.sacc$cellType <- droplevels(sce.sacc$cellType)
## Add cellType.broad
sce.sacc$cellType.Broad <- ss(as.character(sce.sacc$cellType), "\\.", 1)
sce.sacc$cellType.Broad <- factor(sce.sacc$cellType.Broad)
## Match rownames
rownames(sce.sacc) <- rowData(sce.sacc)$ID
table(rownames(sce.sacc) %in% rownames(rse_gene))
# FALSE  TRUE 
# 15021 18517 

## filter for common between sc and bulk data
common_genes <- rownames(rse_gene)[rownames(rse_gene) %in% rownames(sce.sacc)]
sce.sacc <- sce.sacc[common_genes, ]

# Remove 0 genes across all nuclei
sce.sacc <- sce.sacc[!rowSums(assay(sce.sacc, "counts"))==0, ] 
dim(sce.sacc)
# [1] 17785  7004

#### Add bp length for RPKM later #### 
rd <- rowData(sce.sacc)
## Import gnomic features
txdb <- makeTxDbFromGFF("/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/GRCh38-3.0.0_premrna/genes/genes.gtf")
g <- genes(txdb)
e <- exonsBy(txdb, by = "gene")

if (!all(names(e) %in% names(g))) {
  warning("Dropping exons with gene ids not present in the gene list")
  e <- e[names(e) %in% names(g)]
}
e2 <- disjoin(e)
g$bp_length <- sum(width(e2))
summary(g$bp_length)

g_sacc <- g[rownames(sce.sacc),]
summary(g_sacc$bp_length)
table(g_sacc$bp_length > 205012)

rowRanges(sce.sacc) <- g_sacc
all(rownames(rd) == rownames(sce.sacc))
rowData(sce.sacc)$Symbol <- rd$Symbol

## Save filtered sce object
save(sce.sacc, file = here("deconvolution","data","sce.sacc_filtered.Rdata"))
