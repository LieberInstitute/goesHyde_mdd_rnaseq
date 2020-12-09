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

## Amyg Data
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose = TRUE)
sce.amy$uniqueID <- paste0(sce.amy$donor, "_", sce.amy$Barcode)
colnames(sce.amy) <- sce.amy$uniqueID

sce.amy <- sce.amy[,sce.amy$cellType != "Ambig.lowNtrxts",]
sce.amy$cellType <- droplevels(sce.amy$cellType)
## Add cellType.broad
sce.amy$cellType.Broad <- ss(as.character(sce.amy$cellType), "\\.", 1)
sce.amy$cellType.Broad <- factor(sce.amy$cellType.Broad)
## Match rownames
rownames(sce.amy) <- rowData(sce.amy)$ID
table(rownames(sce.amy) %in% rownames(rse_gene))
# FALSE  TRUE 
# 15021 18517 

## filter for common between sc and bulk data
common_genes <- rownames(rse_gene)[rownames(rse_gene) %in% rownames(sce.amy)]
sce.amy <- sce.amy[common_genes, ]

# Remove 0 genes across all nuclei
sce.amy <- sce.amy[!rowSums(assay(sce.amy, "counts"))==0, ] 
dim(sce.amy)
# [1] 17792  6582

#### Add bp length for RPKM later #### 
rd.sacc <- rowData(sce.sacc)
rd.amy <- rowData(sce.amy)
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

## Add data to sacc
g_sacc <- g[rownames(sce.sacc),]
summary(g_sacc$bp_length)
table(g_sacc$bp_length > 205012)

rowRanges(sce.sacc) <- g_sacc
all(rownames(rd.sacc) == rownames(sce.sacc))
rowData(sce.sacc)$Symbol <- rd.sacc$Symbol

## Add data to amy data
g_amy <- g[rownames(sce.amy),]
summary(g_amy$bp_length)
table(g_amy$bp_length > 205012)

rowRanges(sce.amy) <- g_amy
all(rownames(rd.amy) == rownames(sce.amy))
rowData(sce.amy)$Symbol <- rd.amy$Symbol

## Save filtered sce object
save(sce.sacc, file = here("deconvolution","data","sce.sacc_filtered.Rdata"))
save(sce.amy, file = here("deconvolution","data","sce.amy_filtered.Rdata"))


