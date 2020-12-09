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
library(patchwork)

source(here("deconvolution","get_mean_ratio.R"))
source(here("deconvolution","findMarkers_1vAll.R"))
#### Load sce data ####
load(here("deconvolution","data","sce.sacc_filtered.Rdata"), verbose = TRUE)
load(here("deconvolution","data","sce.amyg_filtered.Rdata"), verbose = TRUE)
#### Find mean ratio for all genes ####
ct <- list(broad = "cellType.Broad", specific = "cellType")
sce <- list(sacc = sce.sacc, amyg = sce.amyg)

sceXct <- cross2(sce,ct)
names(sceXct) <- map(cross2(names(sce),names(ct)), ~paste0(.x[[1]],"_",.x[[2]]))
# test <- get_mean_ratio(sce.sacc[1:100,], "cellType.Broad")
mean_ratio <- map(sceXct, ~get_mean_ratio(sce = .x[[1]], cellType_col = .x[[2]]))
map_int(mean_ratio, nrow)
# sacc_broad    amyg_broad sacc_specific amyg_specific 
# 18581         20330         44783         52444 

#### Run findMarkers ####
markers.t.1vAll <- map(sceXct, ~findMarkers_1vAll(sce = .x[[1]], cellType_col = .x[[2]]))
map_int(markers.t.1vAll, nrow)
# sacc_broad    amyg_broad sacc_specific amyg_specific 
# 106710        106752        177850        213504 

#### join stats ####
marker_stats <- map2(mean_ratio, markers.t.1vAll, ~left_join(.x, .y, by = c("gene", "cellType.target","Symbol")) %>%
                       mutate(anno = paste0(" ",anno_ratio,"\n",anno_logFC))
                     )

#### Save tables ####
save(marker_stats, file = here("deconvolution","data","marker_stats.Rdata"))

top_marker_tables <- map2(marker_stats,names(marker_stats), ~filter(.x,rank_ratio<= 5) %>% 
                            select(cellType.target,rank_ratio, gene, Symbol, mean_logcount.target, cellType.nextHighest = cellType, mean_logcount.nextHighest = mean_logcount,ratio, std.logFC) %>%
                            arrange(cellType.target) %>%
                            mutate(genecards_url = paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", Symbol),
                                   region = ss(.y,"_"))
)

walk2(top_marker_tables,names(marker_stats), ~write_csv(.x, paste0("data/top5_markers_sACC-",.y,".csv")))  


#### Create Ratio plots #### 
# load(here("deconvolution","data","marker_stats_sacc.Rdata"), verbose = TRUE)
load("data/cell_colors.Rdata", verbose = TRUE)

ratio_plots <- map2(marker_stats, names(marker_stats), 
                   ~ggplot(.x, aes(x=ratio, y=std.logFC, color = cellType, shape = rank_ratio <= 5))+
                     geom_point() +
                     facet_wrap(~cellType.target, scales = "free_x") +
                     labs(x = "target + 0.01/highest non-target + 0.01",
                          title = paste(.y,"cell types"))+ 
                     scale_color_manual(values = cell_colors) +
                     theme_bw()
)
walk2(ratio_plots, names(ratio_plots), ~ggsave(.x, filename = paste0("plots/expr/ratio_vs_stdFC-",.y,".png"), width = 10))


#### Plot expression ####
## rank
top_n <- 5
plots_ratio <- map2(marker_stats[c(1,3)], names(marker_stats)[c(1,3)], function(x,y){
  cells <- unique(x$cellType.target)
  cellType_col <- ct[[ss(y,"_",2)]]
  plots_ratio <- list()
  for(i in cells){
    marker_stat_cell <- x %>% filter(cellType.target == i, rank_ratio <= top_n) %>%
      mutate(Feature = feature_ratio)
    temp_sce <- sce.sacc[marker_stat_cell$gene,]
    rownames(temp_sce) <- marker_stat_cell$feature_ratio
    
    pe <- plotExpression(temp_sce, exprs_values = "logcounts", features = marker_stat_cell$Feature,
                     x=cellType_col, colour_by=cellType_col, point_alpha=0.5, point_size=.2, ncol=5,
                     add_legend=F) +
        scale_color_manual(values = cell_colors)+
        stat_summary(fun = mean, fun.min = mean, fun.max = mean,
                     geom = "crossbar", width = 0.3) +
        geom_text(data = marker_stat_cell, aes(x = -Inf, y = Inf, label = anno),
                  vjust = "inward", hjust = "inward",size = 2.5)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle(label=paste(y,i, "top", top_n ,"markers: Ranked by mean ratios"))
    
    png_name <- paste0("plots/expr/",y,"_",i,"_top",top_n,"_expression_ratio.png")
    ggsave(plot = pe, filename = png_name, width = 9, height = 3)
    plots <- c(plots, pe)
  }
})

walk2(marker_stats[c(1,3)], names(marker_stats)[c(1,3)], function(x,y){
  cells <- unique(x$cellType.target)
  cellType_col <- ct[[ss(y,"_",2)]]
  # pdf(paste0("plots/expr/sACC-",y,"_top",top_n,"_expression.pdf"))
  # par(mfrow = c(3, 1))
  # plots_ratio <- list()
  for(i in cells){
    marker_stat_cell <- x %>% filter(cellType.target == i, rank_marker <= top_n) %>%
      mutate(Feature = feature_marker)
    temp_sce <- sce.sacc[marker_stat_cell$gene,]
    rownames(temp_sce) <- marker_stat_cell$Feature
    
    pe <- plotExpression(temp_sce, exprs_values = "logcounts", features = marker_stat_cell$Feature,
                         x=cellType_col, colour_by=cellType_col, point_alpha=0.5, point_size=.2, ncol=5,
                         add_legend=F) +
      scale_color_manual(values = cell_colors)+
      stat_summary(fun = mean, fun.min = mean, fun.max = mean,
                   geom = "crossbar", width = 0.3) +
      geom_text(data = marker_stat_cell, aes(x = -Inf, y = Inf, label = anno),
                vjust = "inward", hjust = "inward",size = 2.5)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ggtitle(label=paste(y,i, "top", top_n ,"markers: Ranked by std.logFC"))
    
    png_name <- paste0("plots/expr/",y,"_",i,"_top",top_n,"_expression_marker.png")
    ggsave(plot = pe, filename = png_name, width = 9, height = 3)
    # plots <- c(plots, pe)
    # message(png_name, " COMPLETE!")
  }
  # names(plots) <- cells
  # return(plots)
  # dev.off()
})


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


