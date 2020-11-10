
library(SummarizedExperiment)
library(SingleCellExperiment)
library(jaffelab)
library(here)
library(pheatmap)
library(scuttle)
library(tidyverse)
library(reshape2)
library(here)

source(here("main_colors.R"))
load("cell_colors.Rdata", verbose = TRUE)

## Load rse_gene data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID
rse_gene_amyg <- rse_gene[,rse_gene$BrainRegion == "Amygdala"]
rse_gene_sacc <- rse_gene[,rse_gene$BrainRegion == "sACC"]

#### Top 40 data ####
top40_amyg <- read.csv("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/tables/top40genesLists_Amyg-n2_cellType.split_SN-LEVEL-tests_May2020.csv")
top40_sacc <- read.csv("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/tables/top40genesLists_sACC-n2_cellType_SN-LEVEL-tests_May2020.csv")

#### sACC data ####
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose = TRUE)
## Exclude Ambig.lowNtrxts
sce.sacc <- sce.sacc[,sce.sacc$cellType != "Ambig.lowNtrxts",]
## Add cellType.broad
sce.sacc$cellType.Broad <- ss(as.character(sce.sacc$cellType), "\\.", 1)
## Match rownames
rownames(sce.sacc) <- rowData(sce.sacc)$ID

## Build gene annotation table
top40_anno_sacc <- top40_sacc %>%
  select(grep("1vAll", colnames(top40_sacc))) %>%
  rownames_to_column("n") %>%
  melt(id.vars = "n") %>%
  separate(variable,"_",into = c("cellType",NA)) %>%
  group_by(value) %>% 
  arrange(cellType)%>%
  summarize(n = length(unique(cellType)),
            cell_top40 = ifelse(n > 1, 
                                paste(c("Multi",unique(gsub(".[0-9]", "", cellType))), collapse = "-"),
                                cellType)) %>%
  ungroup() %>%
  dplyr::rename(Symbol = value) %>%
  filter(Symbol %in% rowData(rse_gene)$Symbol,
         Symbol %in% rowData(sce.sacc)$Symbol) %>%
  arrange(cell_top40)

dim(top40_anno_sacc)
# [1] 333   3

top40_anno_sacc$ensemblID <- rowData(rse_gene)$ensemblID[match(top40_anno_sacc$Symbol, rowData(rse_gene)$Symbol)]
top40_anno_sacc %>% filter(n>1) %>% count(cell_top40)
# cell_top40          n
# # A tibble: 3 x 2
# cell_top40            n
# <chr>             <int>
# 1 Multi-Excit           7
# 2 Multi-Excit-Inhib     1
# 3 Multi-Inhib           4

## sacc
top40_anno_sacc <- top40_anno_sacc %>% 
  filter(ensemblID %in% rownames(sce.sacc),
         ensemblID %in% rownames(rse_gene_sacc))

top40all_ensm_sacc <- top40_anno_sacc$ensemblID
length(top40all_ensm_sacc)
# [1] 332

# filter expression data
rse_gene_sacc <- rse_gene_sacc[top40all_ensm_sacc,]
sce.sacc <- sce.sacc[top40all_ensm_sacc,]


#### pseudobulk single cell data ####

pb_sacc <- aggregateAcrossCells(sce.sacc, 
                                id = colData(sce.sacc)[,c("donor","cellType")])
colnames(colData(pb_sacc))[[20]] <- "cellType_num"
colnames(pb_sacc) <- paste0(pb_sacc$cellType, "_", pb_sacc$donor)

dim(pb_sacc)
# [1] 332  20

#rearrange pb_sacc to match top40 anno 
pb_sacc <- pb_sacc[top40_anno_sacc$ensemblID,]
rse_gene_sacc <- rse_gene_sacc[top40_anno_sacc$ensemblID,]

## Transform counts by log and scale
pb_mat_sacc <- as.matrix(assays(pb_sacc)$counts)
pb_mat_log_sacc <- log(pb_mat_sacc + 1)

bulk_mat_sacc <- as.matrix(assays(rse_gene_sacc)$counts)
bulk_mat_log_sacc <- log(bulk_mat_sacc + 1)

## define annotation tables
anno_sc_sacc <- as.data.frame(colData(pb_sacc)[,c("donor","cellType","cellType.Broad")])
anno_bulk_sacc <- as.data.frame(colData(rse_gene_sacc)[,c("PrimaryDx","Experiment")])
## create gene anno obj
gene_anno <- top40_anno_sacc %>% select(cell_top40) %>% as.data.frame()
rownames(gene_anno) <- top40_anno_sacc$ensemblID


## define color scale
breaks_log <- seq(0, max(c(pb_mat_log_sacc, bulk_mat_log_sacc)), length = 100)
breaks_scale <- seq(min(c(pb_mat_sacc, bulk_mat_sacc)), max(c(pb_mat_sacc, bulk_mat_sacc)), length = 100)

## define donor colrs
donor_colors <- c(Br5161 = "black",
                  Br5212 = "lightgrey")

## create annotation color schemes
my_anno_colors <- list(cellType = cell_colors[names(cell_colors) %in% anno_sc_sacc$cellType],
                       cellType.Broad = cell_colors[names(cell_colors) %in% anno_sc_sacc$cellType.Broad],
                       PrimaryDx = mdd_Dx_colors, 
                       Experiment = mdd_dataset_colors, 
                       cell_top40 = cell_colors[names(cell_colors) %in% gene_anno$cell_top40],
                       donor = donor_colors)
#### Create Heat maps ####
## log plots
## Unsorted genes

png("plots/heatmap_log_unclustered_sc_sacc.png", height = 800, width = 580)
pheatmap(pb_mat_log_sacc,
         show_rownames = FALSE,
         show_colnames = FALSE,
         breaks = breaks_log,
         annotation_row = gene_anno,
         annotation_col = anno_sc_sacc, 
         annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "sACC single cell ref - log(counts + 1)")
dev.off()

png("plots/heatmap_log_unclustered_bulk_sacc.png", height = 800, width = 580)
pheatmap(bulk_mat_log_sacc,
         show_rownames = FALSE,
         show_colnames = FALSE,
         breaks = breaks_log,
         annotation_row = gene_anno,
         annotation_col = anno_bulk_sacc,
         annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "sACC bulk - log(counts + 1)")
dev.off()

## Sort genes by sc clustering
png("plots/heatmap_log_sc_sacc.png", height = 800, width = 580)
sc_log_heatmap <- pheatmap(pb_mat_log_sacc,
                           show_rownames = FALSE,
                           show_colnames = FALSE,
                           breaks = breaks_log,
                           annotation_row = gene_anno,
                           annotation_col = anno_sc_sacc, 
                           annotation_colors = my_anno_colors,
                           main = "sACC single cell ref - log(counts + 1)")
dev.off()

bulk_mat_log_sacc <- bulk_mat_log_sacc[sc_log_heatmap$tree_row$order,]

png("plots/heatmap_log_bulk_sacc.png", height = 800, width = 580)
pheatmap(bulk_mat_log_sacc,
         show_rownames = FALSE,
         show_colnames = FALSE,
         breaks = breaks_log,
         annotation_col = anno_bulk_sacc, 
         annotation_row = gene_anno,
         annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "sACC bulk- log(counts + 1)")
dev.off()


#### Amyg Data ####
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose = TRUE)
## Exclude Ambig.lowNtrxts
sce.amy <- sce.amy[,sce.amy$cellType.split != "Ambig.lowNtrxts",]
## Add cellType.broad
sce.amy$cellType.Broad <- ss(as.character(sce.amy$cellType.split), "\\.", 1)
## Match rownames
rownames(sce.amy) <- rowData(sce.amy)$ID

## Build gene annotation table
top40_anno_amyg <- top40_amyg %>%
  select(grep("1vAll", colnames(top40_amyg))) %>%
  rownames_to_column("n") %>%
  melt(id.vars = "n") %>%
  separate(variable,"_",into = c("cellType",NA)) %>%
  group_by(value) %>% 
  arrange(cellType)%>%
  summarize(n = length(unique(cellType)),
            cell_top40 = ifelse(n > 1, 
                                paste(c("Multi",unique(gsub(".[0-9]", "", cellType))), collapse = "-"),
                                cellType)) %>%
  ungroup() %>%
  dplyr::rename(Symbol = value) %>%
  filter(Symbol %in% rowData(rse_gene)$Symbol,
         Symbol %in% rowData(sce.amy)$Symbol) %>%
  arrange(cell_top40)

dim(top40_anno_amyg)
# [1] 377   3

top40_anno_amyg$ensemblID <- rowData(rse_gene)$ensemblID[match(top40_anno_amyg$Symbol, rowData(rse_gene)$Symbol)]
top40_anno_amyg %>% filter(n>1) %>% count(cell_top40)
# # A tibble: 12 x 2
# A tibble: 3 x 2
# cell_top40            n
# <chr>             <int>
# 1 Multi-Excit           4
# 2 Multi-Excit-Inhib     1
# 3 Multi-Inhib          28

## amyg
top40_anno_amyg <- top40_anno_amyg %>% 
  filter(ensemblID %in% rownames(sce.amy),
         ensemblID %in% rownames(rse_gene_amyg))

top40all_ensm_amyg <- top40_anno_amyg$ensemblID
length(top40all_ensm_amyg)
# [1] 376

# filter expression data
rse_gene_amyg <- rse_gene_amyg[top40all_ensm_amyg,]
sce.amy <- sce.amy[top40all_ensm_amyg,]

#### pseudobulk single cell data ####

pb_amyg <- aggregateAcrossCells(sce.amy, 
                                id = colData(sce.amy)[,c("donor","cellType.split")])
colnames(colData(pb_amyg))[[21]] <- "cellType.split_num"
colnames(pb_amyg) <- paste0(pb_amyg$cellType.split, "_", pb_amyg$donor)

dim(pb_amyg)
# [1] 376  21

#rearrange pb_amyg to match top40 anno 
pb_amyg <- pb_amyg[top40_anno_amyg$ensemblID,]
rse_gene_amyg <- rse_gene_amyg[top40_anno_amyg$ensemblID,]

## Transform counts by log and scale
pb_mat_amyg <- as.matrix(assays(pb_amyg)$counts)
pb_mat_log_amyg <- log(pb_mat_amyg + 1)
pb_mat_amyg <- scale(pb_mat_amyg)

bulk_mat_amyg <- as.matrix(assays(rse_gene_amyg)$counts)
bulk_mat_log_amyg <- log(bulk_mat_amyg + 1)
bulk_mat_amyg <- scale(bulk_mat_amyg)

## define annotation tables
anno_sc_amyg <- as.data.frame(colData(pb_amyg)[,c("donor","cellType.split","cellType.Broad")])
anno_bulk_amyg <- as.data.frame(colData(rse_gene_amyg)[,c("PrimaryDx","Experiment")])
## create gene anno obj
gene_anno <- top40_anno_amyg %>% select(cell_top40) %>% as.data.frame()
rownames(gene_anno) <- top40_anno_amyg$ensemblID

## define color scale
breaks_log <- seq(0, max(c(pb_mat_log_amyg, bulk_mat_log_amyg)), length = 100)
breaks_scale <- seq(min(c(pb_mat_amyg, bulk_mat_amyg)), max(c(pb_mat_amyg, bulk_mat_amyg)), length = 100)

## create annotation color schemes

my_anno_colors <- list(cellType.split = cell_colors[names(cell_colors) %in% anno_sc_amyg$cellType.split],
                       cellType.Broad = cell_colors[names(cell_colors) %in% anno_sc_amyg$cellType.Broad],
                       PrimaryDx = mdd_Dx_colors, 
                       Experiment = mdd_dataset_colors, 
                       cell_top40 = cell_colors[names(cell_colors) %in% gene_anno$cell_top40],
                       donor = donor_colors)

#### Create Heat maps ####
## log plots
## Unsorted genes

png("plots/heatmap_log_unclustered_sc_amyg.png", height = 800, width = 580)
pheatmap(pb_mat_log_amyg,
         show_rownames = FALSE,
         show_colnames = FALSE,
         breaks = breaks_log,
         annotation_row = gene_anno,
         annotation_col = anno_sc_amyg, 
         annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "amyg single cell ref - log(counts + 1)")
dev.off()

png("plots/heatmap_log_unclustered_bulk_amyg.png", height = 800, width = 580)
pheatmap(bulk_mat_log_amyg,
         show_rownames = FALSE,
         show_colnames = FALSE,
         breaks = breaks_log,
         annotation_row = gene_anno,
         annotation_col = anno_bulk_amyg,
         annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "amyg bulk - log(counts + 1)")
dev.off()

## Sort genes by sc clustering
png("plots/heatmap_log_sc_amyg.png", height = 800, width = 580)
sc_log_heatmap <- pheatmap(pb_mat_log_amyg,
                           show_rownames = FALSE,
                           show_colnames = FALSE,
                           breaks = breaks_log,
                           annotation_row = gene_anno,
                           annotation_col = anno_sc_amyg, 
                           annotation_colors = my_anno_colors,
                           main = "amyg single cell ref - log(counts + 1)")
dev.off()

bulk_mat_log_amyg <- bulk_mat_log_amyg[sc_log_heatmap$tree_row$order,]

png("plots/heatmap_log_bulk_amyg.png", height = 800, width = 580)
pheatmap(bulk_mat_log_amyg,
         show_rownames = FALSE,
         show_colnames = FALSE,
         breaks = breaks_log,
         annotation_col = anno_bulk_amyg, 
         annotation_row = gene_anno,
         annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "amyg bulk- log(counts + 1)")
dev.off()

# sgejobs::job_single('deconvo_heatmap', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript deconvo_heatmap.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessionsession_info()

