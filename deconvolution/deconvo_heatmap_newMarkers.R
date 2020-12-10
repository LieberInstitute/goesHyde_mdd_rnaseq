
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
library(RColorBrewer)
library(rlang)

# Load colors
source(here("main_colors.R"))
load(here("deconvolution","data","cell_colors.Rdata"), verbose = TRUE)

## Load rse_gene data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID
rse_gene_amyg <- rse_gene[,rse_gene$BrainRegion == "Amygdala"]
rse_gene_sacc <- rse_gene[,rse_gene$BrainRegion == "sACC"]

#### Marker genes ####
load(here("deconvolution","data","marker_stats.Rdata"), verbose = TRUE)

top_n <- 5
marker_genes <- map(marker_stats, ~.x %>%
                      arrange(-ratio) %>%
                      filter(rank_ratio <= top_n) %>%
                      arrange(cellType.target, rank_ratio) %>% 
                      pull(gene))

map_int(marker_genes, length)
# sacc_broad    amyg_broad sacc_specific amyg_specific 
# 30            30            50            60 

#### sce data ####
load(here("deconvolution","data","sce.sacc_filtered.Rdata"), verbose = TRUE)
load(here("deconvolution","data","sce.amyg_filtered.Rdata"), verbose = TRUE)

ct <- list(broad = "cellType.Broad", specific = "cellType")
sce <- list(sacc = sce.sacc, amyg = sce.amyg)

sceXct <- cross2(sce,ct)
names(sceXct) <- map(cross2(names(sce),names(ct)), ~paste0(.x[[1]],"_",.x[[2]]))

## subset genes 
sce_filter <- map2(sceXct, marker_genes, ~.x[[1]][.y,])

##pseudobulk single cell data
pb <- map2(sce_filter,rep(ct,each = 2),function(sce, c){
    sce$pb <- paste0(sce$donor,"_",sce[[c]])
    clusIndex = splitit(sce$pb)
    pbcounts <- sapply(clusIndex, function(ii){
      rowSums(assays(sce)$counts[ ,ii])
    }
    )
    # Compute LSFs at this level
    sizeFactors.PB.all  <- librarySizeFactors(pbcounts)
    # Normalize with these LSFs
    geneExprs.temp <- t(apply(pbcounts, 1, function(x) {log2(x/sizeFactors.PB.all + 1)}))
    rownames(geneExprs.temp) <- rowData(sce)$Symbol
    return(geneExprs.temp)
  }
)

# ## Compute RPKM
# pb_rpkm_sacc <- map(pb_sacc, ~log(getRPKM(.x, "bp_length") +1 )) 
# ## Use symbol as rownames
# rownames(pb_rpkm_sacc[[1]]) <- rowData(sce.sacc)[rownames(pb_rpkm_sacc[[1]]),]$Symbol
# rownames(pb_rpkm_sacc[[2]]) <- rowData(sce.sacc)[rownames(pb_rpkm_sacc[[2]]),]$Symbol

#### Prep for heatmaps ####
## define annotation tables

pb_anno <- map2(sce_filter,rep(ct,each = 2), ~colData(.x) %>% 
                        as_tibble() %>% 
                        select(donor, !!sym(.y)) %>%
                        mutate(pb = paste0(donor, "_",!!sym(.y))) %>%
                        unique() %>%
                        column_to_rownames(var = "pb")
                   )
## create gene anno obj

marker_anno <- map(marker_stats, ~.x %>%
                     arrange(-ratio) %>%
                     filter(rank_ratio <= top_n) %>%
                     mutate(rank_ratio = as.character(rank_ratio)) %>%
                     column_to_rownames(var = "Symbol") %>%
                     select(cellType.target, rank_ratio) 
)


## define donor colrs
donor_colors <- c(Br5161 = "black",
                  Br5212 = "lightgrey")

rank_colors <- brewer.pal(n= top_n, name = "Greens")
names(rank_colors) <- as.character(1:top_n)

my_anno_colors <- map(pb_anno, ~list(cellType.target = cell_colors[levels(.x[,2])],
                                         cellType.Broad = cell_colors[levels(.x[,2])],
                                         cellType = cell_colors[levels(.x[,2])],
                                         donor = donor_colors,
                                         rank_ratio = rank_colors))

#### Create Heat maps - sACC RPKM ####

for(n in names(pb)){
  png(paste0("plots/heatmap_",n,"_top",top_n,".png"), height = 1000, width = 580)
  pheatmap(pb[[n]],
           show_colnames = FALSE,
           annotation_row = marker_anno[[n]],
           annotation_col = pb_anno[[n]], 
           annotation_colors = my_anno_colors[[n]],
           cluster_rows = TRUE,
           main = paste(n,"single cell- log2(pseudobulk counts/size factor) - top",top_n, "ratio rank")
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
