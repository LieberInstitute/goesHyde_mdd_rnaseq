
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
library(pheatmap)
library(sessioninfo)
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
marker_genes_ratio <- map(marker_stats, ~.x %>%
                      arrange(-ratio) %>%
                      filter(rank_ratio <= top_n) %>%
                      arrange(cellType.target, rank_ratio) %>% 
                      pull(gene))

marker_genes_mark <- map(marker_stats, ~.x %>%
                            arrange(-ratio) %>%
                            filter(rank_marker <= top_n) %>%
                            arrange(cellType.target, rank_ratio) %>% 
                            pull(gene))

map2_int(marker_genes_ratio, marker_genes_mark, ~sum(.x %in% .y))
# sacc_broad    amyg_broad sacc_specific amyg_specific 
# 3             4             8            11 

marker_genes_mad <- map(marker_stats,~summarize(.x,
                                      mean_ratio = mean(ratio),
                                      mad_ratio = mad(ratio),
                                      mad_cutoff = mean_ratio + (5*mad_ratio))%>%
                           right_join(.x, by = "cellType.target") %>%
                           filter(ratio > mad_cutoff) %>%
                           pull(gene))
map_int(marker_genes_mad, length)

#### sce data ####
load(here("deconvolution","data","sce.sacc_filtered.Rdata"), verbose = TRUE)
load(here("deconvolution","data","sce.amyg_filtered.Rdata"), verbose = TRUE)

ct <- list(broad = "cellType.Broad", specific = "cellType")
sce <- list(sacc = sce.sacc, amyg = sce.amyg)

sceXct <- cross2(sce,ct)
names(sceXct) <- map(cross2(names(sce),names(ct)), ~paste0(.x[[1]],"_",.x[[2]]))

## subset genes 
sce_filter_ratio <- map2(sceXct, marker_genes_ratio, ~.x[[1]][.y,])
sce_filter_mark  <- map2(sceXct, marker_genes_mark, ~.x[[1]][.y,])
sce_filter_mad <- map2(sceXct, marker_genes_mad, ~.x[[1]][.y,])
##pseudobulk single cell data
pb_ratio <- map2(sce_filter_ratio,rep(ct,each = 2),function(sce, c){
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

pb_mark <- map2(sce_filter_mark,rep(ct,each = 2),function(sce, c){
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

pb_mad <- map2(sce_filter_mad,rep(ct,each = 2),function(sce, c){
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
  rownames(geneExprs.temp) <- rowData(sce)$gene_id
  return(geneExprs.temp)
}
)

map(pb_mad, dim)

#### Prep for heatmaps ####
## define annotation tables

pb_anno <- map2(sce_filter_ratio,rep(ct,each = 2), ~colData(.x) %>% 
                        as_tibble() %>% 
                        select(donor, !!sym(.y)) %>%
                        mutate(pb = paste0(donor, "_",!!sym(.y))) %>%
                        unique() %>%
                        column_to_rownames(var = "pb")
                   )
## create gene anno obj

marker_anno_ratio <- map(marker_stats, ~.x %>%
                     filter(rank_ratio <= top_n) %>%
                     mutate(rank_ratio = as.character(rank_ratio)) %>%
                     column_to_rownames(var = "Symbol") %>%
                     select(cellType.target, rank_ratio) 
)

marker_anno_mark <- map(marker_stats[1:3], ~.x %>%
                           filter(rank_marker <= top_n) %>%
                           mutate(rank_marker = as.character(rank_marker)) %>%
                           select(cellType.target, rank_marker, Symbol) %>% 
                           column_to_rownames(var = "Symbol")
)

marker_anno_mad <- map(marker_stats,~summarize(.x,
                                               mean_ratio = mean(ratio),
                                               mad_ratio = mad(ratio),
                                               mad_cutoff = mean_ratio + (5*mad_ratio))%>%
                         right_join(.x, by = "cellType.target") %>%
                         filter(ratio > mad_cutoff) %>% 
                         select(gene, cellType.target)%>% 
                         column_to_rownames(var = "gene"))


## Duplicate?!
map(marker_stats, ~.x %>%
      filter(rank_marker <= top_n) %>%
      mutate(rank_marker = as.character(rank_marker)) %>%
      select(cellType.target, rank_ratio, Symbol) %>%
      ungroup() %>% count(Symbol) %>% filter(n>1))

## define donor colrs
donor_colors <- c(Br5161 = "black",
                  Br5212 = "lightgrey")

rank_colors <- brewer.pal(n= top_n, name = "Greens")
names(rank_colors) <- as.character(1:top_n)

my_anno_colors <- map(pb_anno, ~list(cellType.target = cell_colors[levels(.x[,2])],
                                         cellType.Broad = cell_colors[levels(.x[,2])],
                                         cellType = cell_colors[levels(.x[,2])],
                                         donor = donor_colors,
                                         rank_ratio = rank_colors,
                                         rank_marker = rank_colors,
                                     PrimaryDx = mdd_Dx_colors,
                                     Experiment = mdd_dataset_colors,
                                     Sex = c(M = "steelblue4", `F` = "mediumvioletred")
                                     )
                      )

#### Create Heat maps - sACC RPKM ####

for(n in names(pb_ratio)){
  png(paste0("plots/heatmap_",n,"_top",top_n,"_ratio.png"), height = 750, width = 550)
  pheatmap(pb_ratio[[n]],
           show_colnames = FALSE,
           annotation_row = marker_anno_ratio[[n]],
           annotation_col = pb_anno[[n]], 
           annotation_colors = my_anno_colors[[n]],
           cluster_rows = TRUE,
           main = paste(n," sn data:log2(norm pseudobulk counts) - top",top_n, "ratio rank")
  )
  dev.off()
}

for(n in names(pb_mark[1:3])){
  png(paste0("plots/heatmap_",n,"_top",top_n,"_mark.png"), height = 750, width = 550)
  pheatmap(pb_mark[[n]],
           show_colnames = FALSE,
           annotation_row = marker_anno_mark[[n]],
           annotation_col = pb_anno[[n]], 
           annotation_colors = my_anno_colors[[n]],
           cluster_rows = TRUE,
           main = paste(n,"sn data:log2(norm pseudobulk counts) - top",top_n, "1vAll rank")
  )
  dev.off()
}

for(n in names(pb_mad)){
  png(paste0("plots/heatmap_",n,"_ratioMAD.png"), height = 750, width = 550)
  pheatmap(pb_mad[[n]],
           show_colnames = FALSE,
           show_rownames = FALSE,
           annotation_row = marker_anno_mad[[n]],
           annotation_col = pb_anno[[n]], 
           annotation_colors = my_anno_colors[[n]],
           cluster_rows = TRUE,
           main = paste(n,"sn data:log2(norm pseudobulk counts) - MAD5")
  )
  dev.off()
}

#### Velmeshev data ####

load("/dcl01/ajaffe/data/lab/singleCell/velmeshev2019/analysis_MNT/SCE_asd-velmeshev-etal_MNT.rda", verbose = TRUE)
map_int(marker_genes_ratio, ~sum(!.x %in% rownames(sce.asd)))
sce.asd <- sce.asd[,sce.asd$region == "ACC"]

cell_type_guess <- data.frame(cluster = levels(sce.asd$cluster),
                              cellType = c(rep("Astro",2),"Micro",rep("Inhib.1",2),
                                           rep("Inhib.2",2),"Excit.1","Excit.2",
                                           "Excit.3","Excit.3","Micro",rep("Other",3),
                                           "Oligo","OPC"))


pb_asd <- map(marker_genes_ratio[c("sacc_broad","sacc_specific")],function(m){
  common_genes <- m[m %in% rownames(sce.asd)]
  sce <- sce.asd[common_genes]
  message(nrow(sce))
  sce$pb <- paste0(sce$individual,"_",sce$cluster)
  clusIndex = splitit(sce$pb)
  pbcounts <- sapply(clusIndex, function(ii){
    rowSums(exp(assays(sce)$logcounts)[ ,ii])
  }
  )
  # Compute LSFs at this level
  sizeFactors.PB.all  <- librarySizeFactors(pbcounts)
  # Normalize with these LSFs
  geneExprs.temp <- t(apply(pbcounts, 1, function(x) {log2(x/sizeFactors.PB.all + 1)}))
  rownames(geneExprs.temp) <- rowData(sce.sacc)[common_genes,]$Symbol
  return(geneExprs.temp)
}
)


marker_anno_vel <- map(marker_stats[c("sacc_broad","sacc_specific")], ~.x %>%
                           filter(rank_ratio <= top_n) %>%
                           mutate(rank_ratio = as.character(rank_ratio)) %>%
                           column_to_rownames(var = "Symbol") %>%
                           select(cellType.target, rank_ratio) 
)

pb_anno_vel <- map(pb_asd, function(x){
  anno <- tibble(id = colnames(x), a = colnames(x)) %>%
    separate(a, "_", into = c("indi","cluster"))%>%
    left_join(cell_type_guess)%>%
    select(-indi) %>%
    column_to_rownames(var = "id")
  return(anno)
})

my_anno_colors_vl <- map2(marker_anno_vel,pb_anno_vel,
                          ~list(cellType.target = c(cell_colors,Other = "grey")[unique(.x[,1])],
                                cellType = c(cell_colors,Other = "grey")[unique(.y[,2])],
                                rank_ratio = rank_colors)
                          
)

breaksList = seq(3, 12, by = 0.08)
for(n in names(pb_asd)){
  png(paste0("plots/heatmap-velmeshev_",n,"_top",top_n,"_ratio.png"), height = 750, width = 550)
  pheatmap(pb_asd[[n]],
           show_colnames = FALSE,
           annotation_row = marker_anno_vel[[n]],
           annotation_col = pb_anno_vel[[n]], 
           annotation_colors = my_anno_colors_vl[[n]],
           cluster_rows = TRUE,
           breaks = breaksList,
           main = paste(n," sn data:log2(norm pseudobulk counts) - top",top_n,"Velmeshev")
  )
  dev.off()
}

#### Bulk Heatmaps ####
bulk_counts <- map2(marker_genes_ratio, rep(list(rse_gene_sacc, rse_gene_amyg),2), function(markers,rse){
  rse_temp <- rse[markers,]
  # rownames(rse_temp) <- ss(rowData(rse_temp)$Symbol,"\\.")
  rownames(rse_temp) <- rowData(sce.sacc)[markers,]$Symbol
  return(assays(rse_temp)$counts)
})

bulk_logcounts <- 

bulk_anno <- as.data.frame(colData(rse_gene)[,c("PrimaryDx",'Experiment','Sex')])

for(n in names(bulk_counts)){
  png(paste0("plots/heatmap-Bulk_",n,"_top",top_n,"_ratio.png"), height = 750, width = 550)
  pheatmap(bulk_counts[[n]],
           show_colnames = FALSE,
           annotation_row = marker_anno_ratio[[n]],
           annotation_col = bulk_anno,
           annotation_colors = my_anno_colors[[n]],
           main = paste(n, "Bulk logcounts"))
  dev.off()
}

# sgejobs::job_single('deconvo_heatmap.R', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript deconvo_heatmap.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
