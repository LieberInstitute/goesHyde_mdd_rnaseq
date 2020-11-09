
library(SummarizedExperiment)
library(SingleCellExperiment)
library(jaffelab)
library(here)
library(pheatmap)
library(scuttle)
library(lattice)
library(tidyverse)
library(reshape2)
library(here)

# source(here("main_colors.R"))
# source("cell_colors.R")
# 
# my_anno_colors <- c(cell_colors,mdd_Dx_colors, mdd_dataset_colors)
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

## Build gene annotation table
top40_anno_sacc <- top40_sacc %>%
  select(grep("1vAll", colnames(top40_sacc))) %>%
  rownames_to_column("n") %>%
  melt(id.vars = "n") %>%
  separate(variable,"_",into = c("cellType",NA)) %>%
  group_by(value) %>% 
  arrange(cellType)%>%
  summarize(n = length(unique(cellType)),
            cell_top40 = paste(cellType, collapse = "_")) %>%
  ungroup() %>%
  dplyr::rename(Symbol = value) %>%
  filter(Symbol %in% rowData(rse_gene)$Symbol &
           Symbol %in% rowData(sce.sacc)$Symbol) %>%
  arrange(cell_top40)

dim(top40_anno_sacc)
# [1] 332   3

top40_anno_sacc$ensemblID <- rowData(rse_gene)$ensemblID[match(top40_anno_sacc$Symbol, rowData(rse_gene)$Symbol)]
## create gene anno obj
gene_anno <- top40_anno_sacc %>% select(cell_top40) %>% as.data.frame()
rownames(gene_anno) <- top40_anno_sacc$ensemblID

top40_anno_sacc %>% filter(n>1) %>% count(cell_top40)
# cell_top40          n
# <chr>           <int>
#   1 Excit.1_Excit.3     2
# 2 Excit.2_Excit.3     2
# 3 Excit.2_Excit.4     1
# 4 Excit.3_Excit.4     2
# 5 Excit.3_Inhib.1     1
# 6 Inhib.1_Inhib.2     4

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

#rearrange pb_sacc to match top40 anno 
pb_sacc <- pb_sacc[rownames(gene_anno),]
rse_gene_sacc <- rse_gene_sacc[row.names(gene_anno),]

## Transform counts by log and scale
pb_mat_sacc <- as.matrix(assays(pb_sacc)$counts)
pb_mat_log_sacc <- log(pb_mat_sacc + 1)
pb_mat_sacc <- scale(pb_mat_sacc)

bulk_mat_sacc <- as.matrix(assays(rse_gene_sacc)$counts)
bulk_mat_log_sacc <- log(bulk_mat_sacc + 1)
bulk_mat_sacc <- scale(bulk_mat_sacc)

## define annotation tables
anno_sc_sacc <- as.data.frame(colData(pb_sacc)[,c("donor", "cellType","cellType.Broad")])
anno_bulk_sacc <- as.data.frame(colData(rse_gene_sacc)[,c("PrimaryDx","Experiment")])

#### Create Heat maps ####
## log plots
## Unsorted genes

pdf("plots/heatmap_log_unclustered_sc_sacc.pdf", height = 10)
pheatmap(pb_mat_log_sacc,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_row = gene_anno,
         annotation_col = anno_sc_sacc, 
         cluster_rows = FALSE,
         main = "sACC single cell ref - log(counts + 1)")
dev.off()

pdf("plots/heatmap_log_unclustered_bulk_sacc.pdf", height = 10)
pheatmap(bulk_mat_log_sacc,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_row = gene_anno,
         annotation_col = anno_bulk_sacc,
         # annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "sACC bulk - log(counts + 1)")
dev.off()

## Sort genes by sc clustering
pdf("plots/heatmap_log_sc_sacc.pdf", height = 10)
sc_log_heatmap <- pheatmap(pb_mat_log_sacc,
                           show_rownames = FALSE,
                           show_colnames = FALSE,
                           annotation_row = gene_anno,
                           annotation_col = anno_sc_sacc, 
                           main = "sACC single cell ref - log(counts + 1)")
dev.off()

bulk_mat_log_sacc <- bulk_mat_log_sacc[sc_log_heatmap$tree_row$order,]

pdf("plots/heatmap_log_bulk_sacc.pdf", height = 10)
pheatmap(bulk_mat_log_sacc,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_col = anno_bulk_sacc, 
         annotation_row = gene_anno,
         cluster_rows = FALSE,
         main = "sACC bulk- log(counts + 1)")
dev.off()

## Scale plots
## Unsorted genes
pdf("plots/heatmap_scale_unclustered_sc_sacc.pdf", height = 10)
pheatmap(pb_mat_sacc,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_row = gene_anno,
         annotation_col = anno_sc_sacc, 
         cluster_rows = FALSE,
         main = "sACC single cell ref")
dev.off()

pdf("plots/heatmap_scale_unclustered_bulk_sacc.pdf", height = 10)
pheatmap(bulk_mat_sacc,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_row = gene_anno,
         annotation_col = anno_bulk_sacc,
         # annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "sACC bulk")
dev.off()

## Sort genes by sc clustering
pdf("plots/heatmap_scale_sc_sacc.pdf", height = 10)
sc_heatmap <- pheatmap(pb_mat_sacc,
                       show_rownames = FALSE,
                       show_colnames = FALSE,
                       annotation_row = gene_anno,
                       annotation_col = anno_sc_sacc, 
                       main = "sACC single cell ref")
dev.off()

bulk_mat_sacc <- bulk_mat_sacc[sc_heatmap$tree_row$order,]

pdf("plots/heatmap_scale_bulk_sacc.pdf", height = 10)
pheatmap(bulk_mat_sacc,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_col = anno_bulk_sacc, 
         annotation_row = gene_anno,
         cluster_rows = FALSE,
         main = "sACC bulk")
dev.off()




