
library(SummarizedExperiment)
library(SingleCellExperiment)
library(jaffelab)
library(here)
library(pheatmap)
library(scuttle)

## Load rse_gene data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID
rse_gene_amyg <- rse_gene[,rse_gene$BrainRegion == "Amygdala"]
rse_gene_sacc <- rse_gene[,rse_gene$BrainRegion == "sACC"]

#### sacc Data ####
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose = TRUE)
## Exclude Ambig.lowNtrxts
sce.sacc <- sce.sacc[,sce.sacc$cellType != "Ambig.lowNtrxts",]
## Add cellType.broad
sce.sacc$cellType.Broad <- ss(as.character(sce.sacc$cellType), "\\.", 1)
## Match rownames
rownames(sce.sacc) <- rowData(sce.sacc)$ID


#### Amyg Data ####
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose = TRUE)
## Exclude Ambig.lowNtrxts
sce.amy <- sce.amy[,sce.amy$cellType.split != "Ambig.lowNtrxts",]
## Add cellType.broad
sce.amy$cellType.Broad <- ss(as.character(sce.amy$cellType.split), "\\.", 1)
## Match rownames
rownames(sce.amy) <- rowData(sce.amy)$ID

#### Top 40 data ####
top40_amyg <- read.csv("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/tables/top40genesLists_Amyg-n2_cellType.split_SN-LEVEL-tests_May2020.csv")
top40_sacc <- read.csv("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/tables/top40genesLists_sACC-n2_cellType_SN-LEVEL-tests_May2020.csv")

## sacc
top40all_sacc <- unique(unlist(top40_sacc[,grepl("1vAll", colnames(top40_sacc))]))
length(top40all_sacc)
# [1] 388
top40all_ensm_sacc <-  rowData(rse_gene)$ensemblID[rowData(rse_gene)$Symbol %in% top40all_sacc]
top40all_ensm_sacc <- top40all_ensm_sacc[top40all_ensm_sacc %in% rownames(sce.sacc)]
length(top40all_ensm_sacc)
# [1] 332

# filter expression data
rse_gene_sacc <- rse_gene_sacc[top40all_ensm_sacc,]
sce.sacc <- sce.sacc[top40all_ensm_sacc,]

##Amyg
top40all_amyg <- unique(unlist(top40_amyg[,grepl("1vAll", colnames(top40_amyg))]))
length(top40all_amyg)
# [1] 441
top40all_ensm_amyg <-  rowData(rse_gene)$ensemblID[rowData(rse_gene)$Symbol %in% top40all_amyg]
top40all_ensm_amyg <- top40all_ensm_amyg[top40all_ensm_amyg %in% rownames(sce.amy)]
length(top40all_ensm_amyg)
# [1] 376  
# filter expression data
rse_gene_amyg <- rse_gene_amyg[top40all_ensm_amyg,]
sce.amy <- sce.amy[top40all_ensm_amyg,]

#### pseudobulk single cell data ####

pb_sacc <- aggregateAcrossCells(sce.sacc, 
                                id = colData(sce.sacc)[,c("donor","cellType")])
