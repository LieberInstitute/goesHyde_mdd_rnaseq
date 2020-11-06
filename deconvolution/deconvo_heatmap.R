
library(SummarizedExperiment)
library(SingleCellExperiment)
library(jaffelab)
library(here)
library(pheatmap)
library(scuttle)
library(lattice)

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
colnames(colData(pb_sacc))[[20]] <- "cellType_num"
colnames(pb_sacc) <- paste0(pb_sacc$cellType, "_", pb_sacc$donor)

dim(pb_sacc)
# [1] 332  20

pb_mat_sacc <- as.matrix(assays(pb_sacc)$counts)
pb_mat_sacc <- scale(pb_mat_sacc)
pb_mat_sacc <- pb_mat_sacc[,order(colnames(pb_mat_sacc))]

bulk_mat_sacc <- as.matrix(assays(rse_gene_sacc)$counts)
bulk_mat_sacc <- scale(bulk_mat_sacc)

anno_sc_sacc <- as.data.frame(colData(pb_sacc)[,c("donor", "cellType","cellType.Broad")])
anno_bulk_sacc <- as.data.frame(colData(rse_gene_sacc)[,c("PrimaryDx","Experiment")])


pdf("plots/heatmap_sc_sacc.pdf")
sc_heatmap <- pheatmap(pb_mat_sacc,
                       show_rownames = FALSE,
                       show_colnames = FALSE,
                       annotation_col = anno_sc_sacc, 
                       main = "sACC single cell ref")
dev.off()

bulk_mat_sacc <- bulk_mat_sacc[sc_heatmap$tree_row$order,]

pdf("plots/heatmap_bulk_sacc.pdf")
pheatmap(bulk_mat_sacc,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_col = anno_bulk_sacc, 
         cluster_rows = FALSE,
         main = "sACC bulk")
dev.off()


