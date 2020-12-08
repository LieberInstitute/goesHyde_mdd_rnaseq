
library(SummarizedExperiment)
library(SingleCellExperiment)
library(jaffelab)
library(here)
library(pheatmap)
library(scuttle)
library(recount)
library(tidyverse)
library(reshape2)
library(here)
library(purrr)

# Load colors
source(here("main_colors.R"))
load(here("deconvolution","data","cell_colors.Rdata"), verbose = TRUE)

## Load rse_gene data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID
rse_gene_amyg <- rse_gene[,rse_gene$BrainRegion == "Amygdala"]
rse_gene_sacc <- rse_gene[,rse_gene$BrainRegion == "sACC"]

#### Marker genes ####
load(here("deconvolution","data","marker_stats_sacc.Rdata"), verbose = TRUE)

top_n <- 10
marker_genes <- map(marker_stats, ~.x %>%
                      arrange(-ratio) %>%
                      slice(1:top_n) %>%
                      arrange(cellType.target, ratio_rank) %>% 
                      pull(gene))

#### sACC data ####
load(here("deconvolution","data","sce.sacc_filtered.Rdata"), verbose = TRUE)

##pseudobulk single cell data

pb_sacc <- map2(marker_genes, c("cellType.Broad", "cellType"), ~aggregateAcrossCells(sce.sacc[.x], 
                                                                                     id = colData(sce.sacc)[,c("donor",.y)]))
map(pb_sacc, dim)
## fix this
colnames(pb_sacc[[1]]) <- paste0(pb_sacc[[1]]$cellType.Broad, "_", pb_sacc[[1]]$donor)
colnames(pb_sacc[[2]]) <- paste0(pb_sacc[[2]]$cellType.Broad, "_", pb_sacc[[2]]$donor)

## Compute RPKM
pb_rpkm_sacc <- map(pb_sacc, ~log(getRPKM(.x, "bp_length") +1 )) 
## Use symbol as rownames
rownames(pb_rpkm_sacc[[1]]) <- rowData(sce.sacc)[rownames(pb_rpkm_sacc[[1]]),]$Symbol
rownames(pb_rpkm_sacc[[2]]) <- rowData(sce.sacc)[rownames(pb_rpkm_sacc[[2]]),]$Symbol
#### Prep for heatmaps ####
## define annotation tables
anno_sc_sacc <- (map2(pb_sacc,c("cellType.Broad","cellType"), ~colData(.x) %>% 
                        as.data.frame() %>% 
                        select(donor, .y)))

## create gene anno obj
marker_anno <- marker_genes_table %>% 
  filter(ratio_rank <= top_n | marker_rank <= top_n) %>%
  select(-ratio_rank, -marker_rank) %>%
  column_to_rownames(var = "gene")

marker_anno <- map(marker_stats, ~.x %>%
                     arrange(-ratio) %>%
                     slice(1:top_n) %>%
                     column_to_rownames(var = "Symbol") %>%
                     select(cellType.target,  ratio) 
)


## define donor colrs
donor_colors <- c(Br5161 = "black",
                  Br5212 = "lightgrey")

cell_colors = cell_colors[unique(c(marker_anno[[1]]$cellType.target,marker_anno[[2]]$cellType.target))]
my_anno_colors <- list(cellType.target = cell_colors,
                       cellType.Broad = cell_colors,
                       cellType = cell_colors,
                       PrimaryDx = mdd_Dx_colors, 
                       Experiment = mdd_dataset_colors, 
                       donor = donor_colors)

#### Create Heat maps - sACC RPKM ####

for(r in names(pb_rpkm_sacc)){
  png(paste0("plots/heatmap_sacc_",r,"_clustered_sc.png"), height = 800, width = 580)
  pheatmap(pb_rpkm_sacc[[r]],
           show_colnames = FALSE,
           annotation_row = marker_anno[[r]],
           annotation_col = anno_sc_sacc[[r]], 
           annotation_colors = my_anno_colors,
           cluster_rows = TRUE,
           main = paste("sACC single cell- log(RPKM + 1) ratio rank - ",r)
  )
  dev.off()
}

# sgejobs::job_single('deconvo_heatmap_newMarkers.R', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript deconvo_heatmap_newMarkers.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessionsession_info()
