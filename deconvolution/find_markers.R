library(SingleCellExperiment)
library(scater)
library(scran)
library(jaffelab)
library(limma)
library(here)
library(reshape2)
library(stringr)
library(purrr)
library(tidyverse)

source(here("deconvolution","get_mean_ratio.R"))
source(here("deconvolution","findMarkers_1vAll.R"))
#### Load sce data ####
load(here("deconvolution","data","sce.sacc_filtered.Rdata"), verbose = TRUE)

#### Find mean ratio for all genes ####
ct <- list(broad = "cellType.Broad", specific = "cellType")

mean_ratio <- map(ct, ~get_mean_ratio(sce = sce.sacc, cellType_col = .x))
map_int(mean_ratio, nrow)
# broad specific 
# 18581    44783 

#### Run findMarkers ####
markers.t.1vAll <- map(ct, ~findMarkers_1vAll(sce = sce.sacc, cellType_col = .x))
map_int(markers.t.1vAll, nrow)
# broad specific 
# 106710   177850

#### join stats ####
marker_stats <- map2(mean_ratio, markers.t.1vAll, ~left_join(.x, .y, by = c("gene", "cellType.target")) %>%
                       filter(log.FDR < -6))
map_int(marker_stats, nrow)
# broad specific 
# 12889    26554
save(marker_stats, file = here("deconvolution","data","marker_stats_sacc.Rdata"))

map(marker_stats, ~count(.x))
# A tibble: 6 x 2
# Groups:   cellType.target [6]
# cellType.target     n
# <chr>           <int>
#   1 Astro             832
# 2 Excit            5162
# 3 Inhib            4186
# 4 Micro             708
# 5 Oligo             841
# 6 OPC              1160
# 
# $specific
# # A tibble: 10 x 2
# # Groups:   cellType.target [10]
# cellType.target     n
# <chr>           <int>
#   1 Astro             832
# 2 Excit.1          4698
# 3 Excit.2          4603
# 4 Excit.3          3355
# 5 Excit.4          2883
# 6 Inhib.1          3887
# 7 Inhib.2          3587
# 8 Micro             708
# 9 Oligo             841
# 10 OPC              1160

#### Create Ratio plots #### 
load("data/cell_colors.Rdata", verbose = TRUE)


ratio_plots <- map2(marker_stats, names(marker_stats), 
                   ~ggplot(.x, aes(x=ratio, y=std.logFC, color = cellType, shape = ratio_rank <= 10))+
                     geom_point(size = .5) +
                     facet_wrap(~cellType.target, scales = "free_x") +
                     labs(x = "target + 0.01/highest non-target + 0.01",
                          title = paste(.y,"cell types"))+ 
                     scale_color_manual(values = cell_colors)
)

ggsave(ratio_plots[[1]] + theme(legend.position = "None") + ratio_plots[[2]], filename = "plots/expr/ratio_vs_stdFC.png", width = 15)

ratio_plots_log10 <- map2(marker_stats, names(marker_stats), 
                    ~ggplot(.x, aes(x=log10(ratio), y=std.logFC, color = cellType, shape = ratio_rank <= 10))+
                      geom_point(size = .5) +
                      facet_wrap(~cellType.target, scales = "free_x") +
                      labs(x = "target + 0.01/highest non-target + 0.01",
                           title = paste(.y,"cell types"))+ 
                      scale_color_manual(values = cell_colors)
)

ggsave(ratio_plots_log10[[1]] + theme(legend.position = "None") + ratio_plots_log10[[2]], filename = "plots/expr/ratio_vs_stdFC_log10.png", width = 15)


mean_plots <- map2(marker_stats, names(marker_stats), 
                   ~ggplot(.x, aes(x=ratio, y=mean_logcount.target, color = cellType, shape = ratio_rank <= 10))+
                     geom_point(size = .5) +
                     facet_wrap(~cellType.target, scales = "free_x") +
                     labs(title = paste(.y,"cell types"))+ 
                     scale_color_manual(values = cell_colors)
)

ggsave(mean_plots[[1]] + theme(legend.position = "None") + mean_plots[[2]], filename = "plots/expr/ratio_vs_mean.png", width = 15)



#### Plot expression ####
for(i in names(markerList.t.1vAll)){
  marker_stat_cell <- marker_stat %>% filter(cellType.marker == i, ratio_rank <= top_n)
  temp_sce <- sce.sacc[marker_stat_cell$gene,]
  rownames(temp_sce) <- marker_stat_cell$Feature
  
  png(paste0("plots/expr/sACC_t-sn-level_1vALL_topMarkers-REFINED-",i,"_logExprs.png"), height=2850, width=1800)
  print(
    plotExpression(temp_sce, exprs_values = "logcounts", features = marker_stat_cell$Feature,
                   x="cellType.Broad", colour_by="cellType.Broad", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + 
      scale_color_manual(values = cell_colors)+
      stat_summary(fun = median, fun.min = median, fun.max = median,
                   geom = "crossbar", width = 0.3) +
      geom_text(data = marker_stat_cell, aes(x = -Inf, y = Inf, label = anno),
                vjust = "inward", hjust = "inward",
                size = 7)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25), text = element_text(size=20)) +  
      # ggtitle(label=paste(i, "top", top_n ,"markers, refined: single-nucleus-level p.w. t-tests, cluster-vs-all")
      ggtitle(label=paste(i, "top", top_n ,"markers: Ranked by median ratios"))
  )
  dev.off()
}

#### Check Overlaps with old data ####
top40_old_sacc <- read.csv("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/tables/top40genesLists_sACC-n2_cellType_SN-LEVEL-tests_May2020.csv")
top40_old_sacc <- top40_old_sacc[,grepl("1vAll", colnames(top40_old_sacc))]

top40_new_sacc <- read.csv("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/tables/top40genesLists-REFINED_sACC-n2_cellType.split_Nov2020.csv")
top40_new_sacc <- top40_new_sacc[,grepl("1vAll", colnames(top40_new_sacc))]

data.frame(cell_types = ss(colnames(top40_old_sacc),"_"),
           overlap_old = map_int(names(top40_old_sacc), ~sum(top40_old_sacc[[.x]] %in% rowData(sce.sacc)[genes.top.t[[ss(ss(.x,"_"),"\\.")]],]$Symbol)),
           overlap_new = map_int(names(top40_new_sacc), ~sum(top40_new_sacc[[.x]] %in% rowData(sce.sacc)[genes.top.t[[ss(ss(.x,"_"),"\\.")]],]$Symbol)))
#    cell_types overlap_old overlap_new
# 1       Astro          40          40
# 2     Excit.1          13          13
# 3     Excit.2           8           8
# 4     Excit.3           3           4
# 5     Excit.4           0           0
# 6     Inhib.1          20          20
# 7     Inhib.2          15          18
# 8       Micro          38          38
# 9       Oligo          39          39
# 10        OPC          37          37 


