library(SingleCellExperiment)
library(jaffelab)
library(here)

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

## Save filtered sce object
save(sce.sacc, file = here("deconvolution","data","sce.sacc_filtered.Rdata"))
