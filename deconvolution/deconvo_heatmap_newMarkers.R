
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
library(purrr)

source(here("main_colors.R"))
load(here("deconvolution","data","cell_colors.Rdata"), verbose = TRUE)

## Load rse_gene data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID
rse_gene_amyg <- rse_gene[,rse_gene$BrainRegion == "Amygdala"]
rse_gene_sacc <- rse_gene[,rse_gene$BrainRegion == "sACC"]

#### Marker genes ####
broad_markers_sacc <- read_csv(here("deconvolution","data","broad_markers_sacc.csv"))
broad_ratio_markers_sacc <- read_csv(here("deconvolution","data","braod_ratio_markers_sacc.csv."))

broad_markers_sacc <- full_join(broad_markers_sacc, broad_ratio_markers_sacc) 
top_n <- 40

rank_type <- list("marker","ratio")
names(rank_type) <- rank_type
marker_genes <- map(rank_type, ~broad_markers_sacc %>%
                      rename(rank = paste0(.x,"_rank")) %>%
                      filter(rank <= 40) %>%
                      arrange(cellType.marker, rank) %>% 
                      pull(gene))

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

g_sacc <- g[broad_markers_sacc$gene,]
summary(g_sacc$bp_length)
table(g_sacc$bp_length > 205012)

rowRanges(sce.sacc) <- g_sacc

#### pseudobulk single cell data ####
## findMarker ranked genes

pb_sacc <- map(marker_genes, )

pb_sacc <- map(marker_genes, ~aggregateAcrossCells(sce.sacc[.x], 
                                                   id = colData(sce.sacc)[,c("donor","cellType.Broad")]))
colnames(pb_sacc[["marker"]])
## fix this
colnames(pb_sacc[[1]]) <- paste0(pb_sacc[[1]]$cellType.Broad, "_", pb_sacc[[1]]$donor)
colnames(pb_sacc[[2]]) <- paste0(pb_sacc[[2]]$cellType.Broad, "_", pb_sacc[[2]]$donor)

## Compute RPKM
pb_rpkm_sacc <- map(pb_sacc, ~log(getRPKM(.x, "bp_length") +1 )) 

#### Prep for heatmaps ####
## define annotation tables
anno_sc_sacc <- as.data.frame(colData(pb_sacc[[1]])[,c("donor","cellType.Broad")])

## create gene anno obj
marker_anno <- broad_markers_sacc %>% 
  filter(marker_rank <= top_n | ratio_rank <= top_n) %>%
  column_to_rownames(var = "gene") %>%
  select(-ratio_rank, -marker_rank) 

## define donor colrs
donor_colors <- c(Br5161 = "black",
                  Br5212 = "lightgrey")

## create annotation color schemes
# my_anno_colors <- list(cellType = cell_colors[names(cell_colors) %in% anno_sc_sacc$cellType],
#                        cellType.Broad = cell_colors[names(cell_colors) %in% anno_sc_sacc$cellType.Broad],
#                        PrimaryDx = mdd_Dx_colors, 
#                        Experiment = mdd_dataset_colors, 
#                        cell_top40 = cell_colors[names(cell_colors) %in% gene_anno$cell_top40],
#                        donor = donor_colors)
cell_colors <- cell_colors[unique(broad_markers_sacc$cellType.marker)]
my_anno_colors <- list(cellType.marker = cell_colors,
                       cellType.Broad = cell_colors,
                       PrimaryDx = mdd_Dx_colors, 
                       Experiment = mdd_dataset_colors, 
                       donor = donor_colors)

#### Create Heat maps - sACC RPKM ####
for(r in rank_type){
  png(paste0("plots/heatmap_sacc_",r,"-rank_unclustered_sc.png"), height = 800, width = 580)
  pheatmap(pb_rpkm_sacc[[r]],
           show_rownames = FALSE,
           show_colnames = FALSE,
           annotation_row = marker_anno,
           annotation_col = anno_sc_sacc, 
           annotation_colors = my_anno_colors,
           cluster_rows = FALSE,
           main = paste("sACC single cell- log(RPKM + 1) rank",r)
  )
  dev.off()
  png(paste0("plots/heatmap_sacc_",r,"-rank_clustered_sc.png"), height = 800, width = 580)
  pheatmap(pb_rpkm_sacc[[r]],
           show_rownames = FALSE,
           show_colnames = FALSE,
           annotation_row = marker_anno,
           annotation_col = anno_sc_sacc, 
           annotation_colors = my_anno_colors,
           cluster_rows = TRUE,
           main = paste("sACC single cell- log(RPKM + 1) rank",r)
  )
  dev.off()
}


