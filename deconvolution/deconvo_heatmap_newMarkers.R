
library(SummarizedExperiment)
library(SingleCellExperiment)
library(jaffelab)
library(here)
library(pheatmap)
library(scuttle)
library(recount)
library(GenomicFeatures)
library(tidyverse)
library(reshape2)
library(here)
library(ggbeeswarm)

source(here("main_colors.R"))
load(here("deconvolution","data","cell_colors.Rdata"), verbose = TRUE)

## Load rse_gene data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID
rse_gene_amyg <- rse_gene[,rse_gene$BrainRegion == "Amygdala"]
rse_gene_sacc <- rse_gene[,rse_gene$BrainRegion == "sACC"]

#### Marker genes ####
broad_markers_sacc <- read.csv(here("deconvolution","data","broad_markers_sacc.csv"))
dim(broad_markers_sacc)
# [1] 240   3

#### sACC data ####
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose = TRUE)
sce.sacc$uniqueID <- paste0(sce.sacc$donor, "_", sce.sacc$Barcode)
colnames(sce.sacc) <- sce.sacc$uniqueID
## Exclude Ambig.lowNtrxts
sce.sacc <- sce.sacc[,sce.sacc$cellType != "Ambig.lowNtrxts",]
sce.sacc$cellType <- factor(sce.sacc$cellType)
## Add cellType.broad
sce.sacc$cellType.Broad <- ss(as.character(sce.sacc$cellType), "\\.", 1)
## Match rownames
rownames(sce.sacc) <- rowData(sce.sacc)$ID

# filter expression data
rse_gene_sacc <- rse_gene_sacc[broad_markers_sacc$gene,]
sce.sacc <- sce.sacc[broad_markers_sacc$gene,]

#### Import gnomic features ####
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

g_sacc <- g[top40all_ensm_sacc,]
summary(g_sacc$bp_length)
table(g_sacc$bp_length > 205012)

rowRanges(sce.sacc) <- g_sacc

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

## Compute RPKM
assays(sce.sacc)$counts <- as.array(assays(sce.sacc)$counts)
sce_rpkm_sacc <- getRPKM(sce.sacc, "bp_length")
pb_rpkm_sacc <- getRPKM(pb_sacc, "bp_length")
bulk_rpkm_sacc <- getRPKM(rse_gene_sacc, "Length")

#### Plot Mean RPKM ####
cell_types_sacc <- as.tibble(colData(sce.sacc)) %>% 
  select(sample = uniqueID , cellType) %>%
  mutate(cellType = as.character(cellType))

rpkm_long_sacc <- cbind(pb_rpkm_sacc, bulk_rpkm_sacc) %>% 
  melt() %>%
  rename(Gene = Var1, sample = Var2, RPKM = value) 
# %>%
# left_join(cell_types_sacc, by = "sample") %>%
# replace_na(list(cellType = "bulk"))

rpkm_long_sacc$cellType <- factor(rpkm_long_sacc$cellType)
rpkm_long_sacc$cellType <- relevel(rpkm_long_sacc$cellType, "bulk")

mean_rpkm_sacc <- rpkm_long_sacc %>% group_by(sample) %>%
  summarise(mean_rpkm = mean(RPKM),
            sum_rpkm = sum(RPKM)) %>%
  ungroup() %>%
  mutate(id = row_number()) %>%
  separate(sample, "_", into = c("cellType", "donor",NA))  %>%
  mutate(cellType = ifelse(grepl('^R', cellType), "bulk",cellType))

scatter_mean_rpkm <- mean_rpkm_sacc%>%
  ggplot(aes(id, mean_rpkm, color= cellType))+
  geom_point()

ggsave(plot = scatter_mean_rpkm, filename = "plots/meanRPKM_scatter_sacc.png")

boxplot_mean_rpkm <- mean_rpkm_sacc%>%
  ggplot(aes(cellType, mean_rpkm, fill= cellType))+
  geom_boxplot() +
  labs( x= "Data type", y = "Sample Mean RPKM", title = "Mean RPKM Distribution",
        subtitle = "sACC") + 
  scale_fill_manual(values = c(cell_colors, bulk = "grey"))

ggsave(plot = boxplot_mean_rpkm, filename = "plots/meanRPKM_boxplot_sacc.png")


boxplot_sum_rpkm <- mean_rpkm_sacc%>%
  ggplot(aes(cellType, sum_rpkm, fill= cellType))+
  geom_boxplot() +
  labs( x= "Data type", y = "Sum of Sample RPKM", title = "Top40 marker genes per cell type",
        subtitle = "sACC") + 
  scale_fill_manual(values = c(cell_colors, bulk = "grey"))

ggsave(plot = boxplot_sum_rpkm, filename = "plots/sumRPKM_boxplot_sacc.png")


beeswarm_mean_rpkm <- mean_rpkm_sacc%>%
  ggplot(aes(cellType, mean_rpkm, color= cellType))+
  geom_beeswarm() +
  labs( x= "Data type", y = "Sample Mean RPKM", title = "Mean RPKM Distribution",
        subtitle = "sACC") + 
  scale_color_manual(values = c(cell_colors, bulk = "grey"))

ggsave(plot = beeswarm_mean_rpkm, filename = "plots/meanRPKM_beeswarm_sacc.png")

rpkm_by_gene <- rpkm_long_sacc %>%
  group_by(Gene, bulk = cellType == "bulk") %>%
  summarise(sd = sd(RPKM), mean = mean(RPKM)) %>% dim()
left_join(top40_anno_sacc %>% select(cell_top40, Gene = ensemblID)) %>%
  mutate(Data = ifelse(bulk, "Bulk","SC"))


rpkm_gene_scatter <- rpkm_by_gene %>% 
  ggplot(aes(mean, sd, color = cell_top40)) +
  geom_point() +
  scale_color_manual(values = cell_colors)+
  facet_wrap(~Data, scales = "free" ) +
  labs(x = 'Mean RPKM', y = 'SD RPKM')

ggsave(plot = rpkm_gene_scatter, filename = "plots/rpkm_gene_scatter_sacc.png", width = 15)

## Transform counts and rpkm by log
pb_counts_log_sacc <- log(as.matrix(assays(pb_sacc)$counts) + 1)
pb_rpkm_log_sacc <- log(pb_rpkm_sacc + 1)
bulk_counts_log_sacc <- log(as.matrix(assays(rse_gene_sacc)$counts) + 1)
bulk_rpkm_log_sacc <- log(bulk_rpkm_sacc +1)

## define color scale
breaks_count <- seq(0, max(c(pb_counts_log_sacc, bulk_counts_log_sacc)), length = 100)
breaks_rpkm <- seq(0, max(c(pb_rpkm_log_sacc, bulk_rpkm_log_sacc)), length = 100)

## define annotation tables
anno_sc_sacc <- as.data.frame(colData(pb_sacc)[,c("donor","cellType","cellType.Broad")])
anno_bulk_sacc <- as.data.frame(colData(rse_gene_sacc)[,c("PrimaryDx","Experiment")])

## create gene anno obj
gene_anno <- top40_anno_sacc %>% select(cell_top40) %>% as.data.frame()
rownames(gene_anno) <- top40_anno_sacc$ensemblID

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
#### Create Heat maps - sacc counts ####
## Unsorted genes
png("plots/heatmap_counts_sacc_sc_unclustered_nov.png", height = 800, width = 580)
pheatmap(pb_counts_log_sacc,
         show_rownames = FALSE,
         show_colnames = FALSE,
         breaks = breaks_count,
         annotation_row = gene_anno,
         annotation_col = anno_sc_sacc, 
         annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "sACC single cell ref - log(counts + 1)")
dev.off()

png("plots/heatmap_counts_sacc_bulk_unclustered.png", height = 800, width = 580)
pheatmap(bulk_counts_log_sacc,
         show_rownames = FALSE,
         show_colnames = FALSE,
         breaks = breaks_count,
         annotation_row = gene_anno,
         annotation_col = anno_bulk_sacc,
         annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "sACC bulk - log(counts + 1)")
dev.off()

## Sort by 
png("plots/heatmap_counts_sacc_bulk_explore.png", height = 800, width = 580)
bulk_explore_sacc <- pheatmap(bulk_counts_log_sacc,
                              show_rownames = FALSE,
                              show_colnames = FALSE,
                              breaks = breaks_count,
                              annotation_row = gene_anno,
                              annotation_col = anno_bulk_sacc,
                              annotation_colors = my_anno_colors,
                              cutree_cols = 4,
                              main = "sACC bulk - log(counts + 1)")
dev.off()

## Sort genes by sc clustering
png("plots/heatmap_counts_sacc_sc_nov.png", height = 800, width = 580)
sc_log_heatmap <- pheatmap(pb_counts_log_sacc,
                           show_rownames = FALSE,
                           show_colnames = FALSE,
                           breaks = breaks_count,
                           annotation_row = gene_anno,
                           annotation_col = anno_sc_sacc, 
                           annotation_colors = my_anno_colors,
                           main = "sACC single cell ref - log(counts + 1)")
dev.off()

bulk_counts_log_sacc <- bulk_counts_log_sacc[sc_log_heatmap$tree_row$order,]

png("plots/heatmap_counts_sacc_bulk.png", height = 800, width = 580)
pheatmap(bulk_counts_log_sacc,
         show_rownames = FALSE,
         show_colnames = FALSE,
         breaks = breaks_count,
         annotation_col = anno_bulk_sacc, 
         annotation_row = gene_anno,
         annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "sACC bulk- log(counts + 1)")
dev.off()

#### Create Heat maps - sACC RPKM ####
## Unsorted genes
png("plots/heatmap_rpkm_sacc_sc_unclustered.png", height = 800, width = 580)
pheatmap(pb_rpkm_log_sacc,
         show_rownames = FALSE,
         show_colnames = FALSE,
         # #breaks = breaks_rpkm,
         annotation_row = gene_anno,
         annotation_col = anno_sc_sacc, 
         annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "sACC single cell ref - log(RPKM + 1)")
dev.off()

png("plots/heatmap_rpkm_sacc_bulk_unclustered.png", height = 800, width = 580)
pheatmap(bulk_rpkm_log_sacc,
         show_rownames = FALSE,
         show_colnames = FALSE,
         #breaks = breaks_rpkm,
         annotation_row = gene_anno,
         annotation_col = anno_bulk_sacc,
         annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "sACC bulk - log(RPKM + 1)")
dev.off()

## Sort by genes by bulk data
png("plots/heatmap_rpkm_sacc_bulk_explore.png", height = 800, width = 580)
bulk_explore_sacc <- pheatmap(bulk_rpkm_log_sacc,
                              show_rownames = FALSE,
                              show_colnames = FALSE,
                              #breaks = breaks_rpkm,
                              annotation_row = gene_anno,
                              annotation_col = anno_bulk_sacc,
                              annotation_colors = my_anno_colors,
                              cutree_cols = 4,
                              main = "sACC bulk - log(rpkm + 1)")
dev.off()

## Sort genes by sc clustering
png("plots/heatmap_rpkm_sacc_sc.png", height = 800, width = 580)
sc_log_heatmap <- pheatmap(pb_rpkm_log_sacc,
                           show_rownames = FALSE,
                           show_colnames = FALSE,
                           #breaks = breaks_rpkm,
                           annotation_row = gene_anno,
                           annotation_col = anno_sc_sacc, 
                           annotation_colors = my_anno_colors,
                           main = "sACC single cell ref - log(rpkm + 1)")
dev.off()

bulk_rpkm_log_sacc <- bulk_rpkm_log_sacc[sc_log_heatmap$tree_row$order,]

png("plots/heatmap_rpkm_sacc_bulk.png", height = 800, width = 580)
pheatmap(bulk_rpkm_log_sacc,
         show_rownames = FALSE,
         show_colnames = FALSE,
         #breaks = breaks_rpkm,
         annotation_col = anno_bulk_sacc, 
         annotation_row = gene_anno,
         annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "sACC bulk- log(rpkm + 1)")
dev.off()

#### Amyg Data ####
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose = TRUE)
## Exclude Ambig.lowNtrxts
sce.amy <- sce.amy[,sce.amy$cellType.split != "Ambig.lowNtrxts",]
sce.amy$cellType.split <- factor(sce.amy$cellType.split)
## Add cellType.broad
sce.amy$cellType.Broad <- ss(as.character(sce.amy$cellType.split), "\\.", 1)
## Match rownames
rownames(sce.amy) <- rowData(sce.amy)$ID

round(table(sce.amy$cellType.split)/ncol(sce.amy),3)
# Astro Excit.1 Excit.2 Excit.3 Inhib.1 Inhib.2 Inhib.3 Inhib.4 Inhib.5   Micro   Oligo     OPC 
# 0.129   0.051   0.006   0.008   0.026   0.017   0.005   0.004   0.015   0.116   0.528   0.095 
round(table(sce.amy$cellType.Broad)/ncol(sce.amy),3)
# Astro Excit Inhib Micro Oligo   OPC 
# 0.129 0.065 0.066 0.116 0.528 0.095 
table(sce.amy$cellType.Broad, sce.amy$donor)

#       Br5161 Br5212
# Astro    489    363
# Excit    141    288
# Inhib    169    268
# Micro    425    339
# Oligo   1697   1776
# OPC      335    292

## Build gene annotation table
top40_anno_amyg <- top40_new_amyg %>%
  dplyr::select(grep("1vAll", colnames(top40_amyg))) %>%
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

## Add rowRanges
g_amyg <- g[top40all_ensm_amyg,]
summary(g_amyg$bp_length)
table(g_amyg$bp_length > 205012)

rowRanges(sce.amy) <- g_amyg

#### pseudobulk single cell data ####

pb_amyg <- aggregateAcrossCells(sce.amy, 
                                id = colData(sce.amy)[,c("donor","cellType.split")])
colnames(colData(pb_amyg))[[21]] <- "cellType.split_num"
colnames(pb_amyg) <- paste0(pb_amyg$cellType.split, "_", pb_amyg$donor)

dim(pb_amyg)
# [1] 376  21

## Compute RPKM
pb_rpkm_amyg <- getRPKM(pb_amyg, "bp_length")
bulk_rpkm_amyg <- getRPKM(rse_gene_amyg, "Length")

## Transform counts by log and scale
pb_counts_log_amyg <- log(as.matrix(assays(pb_amyg)$counts) + 1)
pb_rpkm_log_amyg <- log(pb_rpkm_amyg + 1)
bulk_counts_log_amyg <- log(as.matrix(assays(rse_gene_amyg)$counts) + 1)
bulk_rpkm_log_amyg <- log(bulk_rpkm_amyg + 1)

## define annotation tables
anno_sc_amyg <- as.data.frame(colData(pb_amyg)[,c("donor","cellType.split","cellType.Broad")])
anno_bulk_amyg <- as.data.frame(colData(rse_gene_amyg)[,c("PrimaryDx","Experiment")])
## create gene anno obj
gene_anno <- top40_anno_amyg %>% select(cell_top40) %>% as.data.frame()
rownames(gene_anno) <- top40_anno_amyg$ensemblID

## define color scale
breaks_count <- seq(0, max(c(pb_counts_log_amyg, bulk_counts_log_amyg)), length = 100)
breaks_rpkm <- seq(0, max(c(pb_rpkm_log_amyg, bulk_rpkm_log_amyg)), length = 100)

## create annotation color schemes

my_anno_colors <- list(cellType.split = cell_colors[names(cell_colors) %in% anno_sc_amyg$cellType.split],
                       cellType.Broad = cell_colors[names(cell_colors) %in% anno_sc_amyg$cellType.Broad],
                       PrimaryDx = mdd_Dx_colors, 
                       Experiment = mdd_dataset_colors, 
                       cell_top40 = cell_colors[names(cell_colors) %in% gene_anno$cell_top40],
                       donor = donor_colors)

#### Create Heat maps - amyg counts ####
## log plots
## Unsorted genes

png("plots/heatmap_counts_amyg_sc_unclustered.png", height = 800, width = 580)
pheatmap(pb_counts_log_amyg,
         show_rownames = FALSE,
         show_colnames = FALSE,
         breaks = breaks_count,
         annotation_row = gene_anno,
         annotation_col = anno_sc_amyg, 
         annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "Amyg single cell ref - log(counts + 1)")
dev.off()

png("plots/heatmap_counts_amyg_bulk_unclustered.png", height = 800, width = 580)
pheatmap(bulk_counts_log_amyg,
         show_rownames = FALSE,
         show_colnames = FALSE,
         breaks = breaks_count,
         annotation_row = gene_anno,
         annotation_col = anno_bulk_amyg,
         annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "Amyg bulk - log(counts + 1)")
dev.off()

png("plots/heatmap_counts_amyg_bulk_explore.png", height = 800, width = 580)
bulk_explore <- pheatmap(bulk_counts_log_amyg,
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         breaks = breaks_count,
                         annotation_row = gene_anno,
                         annotation_col = anno_bulk_amyg,
                         annotation_colors = my_anno_colors,
                         cutree_cols = 4,
                         main = "Amyg bulk - log(counts + 1)")
dev.off()

## Sort genes by sc clustering
png("plots/heatmap_counts_amyg_sc.png", height = 800, width = 580)
sc_log_heatmap <- pheatmap(pb_counts_log_amyg,
                           show_rownames = FALSE,
                           show_colnames = FALSE,
                           breaks = breaks_count,
                           annotation_row = gene_anno,
                           annotation_col = anno_sc_amyg, 
                           annotation_colors = my_anno_colors,
                           main = "Amyg single cell ref - log(counts + 1)")
dev.off()

bulk_counts_log_amyg <- bulk_counts_log_amyg[sc_log_heatmap$tree_row$order,]

png("plots/heatmap_counts_amyg_bulk.png", height = 800, width = 580)
pheatmap(bulk_counts_log_amyg,
         show_rownames = FALSE,
         show_colnames = FALSE,
         breaks = breaks_count,
         annotation_col = anno_bulk_amyg, 
         annotation_row = gene_anno,
         annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "amyg bulk- log(counts + 1)")
dev.off()

#### Create Heat maps - amyg rpkm ####
## log plots
## Unsorted genes

png("plots/heatmap_rpkm_amyg_sc_unclustered.png", height = 800, width = 580)
pheatmap(pb_rpkm_log_amyg,
         show_rownames = FALSE,
         show_colnames = FALSE,
         #breaks = breaks_rpkm,
         annotation_row = gene_anno,
         annotation_col = anno_sc_amyg, 
         annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "Amyg single cell ref - log(RPKM + 1)")
dev.off()

png("plots/heatmap_rpkm_amyg_bulk_unclustered.png", height = 800, width = 580)
pheatmap(bulk_rpkm_log_amyg,
         show_rownames = FALSE,
         show_colnames = FALSE,
         #breaks = breaks_rpkm,
         annotation_row = gene_anno,
         annotation_col = anno_bulk_amyg,
         annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "Amyg bulk - log(RPKM + 1)")
dev.off()

png("plots/heatmap_rpkm_amyg_bulk_explore.png", height = 800, width = 580)
bulk_explore <- pheatmap(bulk_rpkm_log_amyg,
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         #breaks = breaks_rpkm,
                         annotation_row = gene_anno,
                         annotation_col = anno_bulk_amyg,
                         annotation_colors = my_anno_colors,
                         cutree_cols = 4,
                         main = "Amyg bulk - log(RPKM + 1)")
dev.off()

## Sort genes by sc clustering
png("plots/heatmap_rpkm_amyg_sc.png", height = 800, width = 580)
sc_log_heatmap <- pheatmap(pb_rpkm_log_amyg,
                           show_rownames = FALSE,
                           show_colnames = FALSE,
                           #breaks = breaks_rpkm,
                           annotation_row = gene_anno,
                           annotation_col = anno_sc_amyg, 
                           annotation_colors = my_anno_colors,
                           main = "Amyg single cell ref - log(RPKM + 1)")
dev.off()

bulk_rpkm_log_amyg <- bulk_rpkm_log_amyg[sc_log_heatmap$tree_row$order,]

png("plots/heatmap_rpkm_amyg_bulk.png", height = 800, width = 580)
pheatmap(bulk_rpkm_log_amyg,
         show_rownames = FALSE,
         show_colnames = FALSE,
         #breaks = breaks_rpkm,
         annotation_col = anno_bulk_amyg, 
         annotation_row = gene_anno,
         annotation_colors = my_anno_colors,
         cluster_rows = FALSE,
         main = "amyg bulk- log(rpkm + 1)")
dev.off()

# sgejobs::job_single('deconvo_heatmap', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript deconvo_heatmap.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessionsession_info()

