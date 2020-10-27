library(RColorBrewer)
library(ggrepel)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(jaffelab)
library(here)
library(MuSiC)
library(Biobase)
library(xbioc)

load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda", verbose = TRUE)
dim(sce.hpc)
# [1] 33538 10444

## drop the "Ambig.lowNtrxts" from sce.hpc$cellType.split
table(sce.hpc$cellType.split ==  "Ambig.lowNtrxts")
# FALSE  TRUE 
# 10343   101
sce.hpc <- sce.hpc[,sce.hpc$cellType.split !=  "Ambig.lowNtrxts"]

## Add cellType.broad
table(sce.hpc$cellType.split)
sce.hpc$cellType.Broad <- ss(as.character(sce.hpc$cellType.split), "\\.", 1)
table(sce.hpc$cellType.Broad)
# Astro Excit Inhib Micro Oligo   OPC Tcell 
# 1343   602   370  1253  5885   864    26


## Load rse_gene data
load(here("data", "rse_gene_GoesZandi.rda"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$Symbol
## Create input Expression sets for MuSiC
rse_gene_sACC <- rse_gene[!duplicated(rownames(rse_gene)),rse_gene$BrainRegion == "sACC"]
es_gene_sACC <- ExpressionSet(assayData = assays(rse_gene_sACC)$counts)

#create unique colnames
sce.hpc$uniqueID <- paste0(sce.hpc$donor, "_", sce.hpc$Barcode)
colnames(sce.hpc) <- sce.hpc$uniqueID
rownames(sce.hpc) <- rowData(sce.hpc)$Symbol
sce.hpc <- sce.hpc[!duplicated(rownames(sce.hpc)),]
# create pheono Data
pd_sce_sacc <- as.data.frame(colData(sce.hpc)[c("cellType.Broad", "uniqueID","sizeFactor","sizeFactor")])
# create 
es_sc_sACC <- ExpressionSet(assayData = as.matrix(assays(sce.hpc)$counts),
                            phenoData=AnnotatedDataFrame(pd_sce_sacc))


table(rownames(es_sc_sACC) %in% rownames(es_gene_sACC))
# FALSE  TRUE 
# 11120 22394 

est_prop_sacc = music_prop(bulk.eset = es_gene_sACC, 
                           sc.eset = es_sc_sACC, 
                           clusters = 'cellType.Broad',
                           samples = 'uniqueID')
save(est_prop_sacc, file = "prop_sacc.Rdata")
