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
library(sessioninfo)

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
# load(here("deconvolution","data","marker_stats.Rdata"), verbose = TRUE)

## save top 5 tables
top_n <- 5
top_marker_tables <- map2(marker_stats,names(marker_stats), ~filter(.x,rank_ratio<= top_n) %>% 
                            select(cellType.target,rank_ratio, gene, Symbol, mean_logcount.target, cellType.nextHighest = cellType, mean_logcount.nextHighest = mean_logcount,ratio, std.logFC) %>%
                            arrange(cellType.target) %>%
                            mutate(genecards_url = paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", Symbol),
                                   region = ss(.y,"_"))
)
## save 
walk2(top_marker_tables,names(marker_stats), ~write_csv(.x, paste0("data/top",top_n,"_markers_",.y,".csv")))  

#### MAD analysis ####
map(marker_stats, ~filter(.x, log.p.value < log(0.05)) %>% summarize(n = n(), min_fc =min(std.logFC)))

map(marker_stats, ~summarize(.x, genes = n(), 
                             n_sig1vAll = sum(log.p.value < log(0.05)),
                             mean_ratio = mean(ratio),
                             ratio_over1 = sum(ratio > 1)))


map(marker_stats, ~summarize(.x,
                             mean_ratio = mean(ratio),
                             mad_ratio = mad(ratio),
                             mad3_cutoff = mean_ratio + (3*mad_ratio),
                             n_markers = sum(ratio > mad3_cutoff)))

## 5 MADs are is the hightest cuttoff to keep some markers for 
map(marker_stats, ~summarize(.x,
                             mean_ratio = mean(ratio),
                             mad_ratio = mad(ratio),
                             mad_cutoff = mean_ratio + (5*mad_ratio),
                             n_markers = sum(ratio > mad_cutoff)))

map(marker_stats, ~filter(.x, log.p.value < log(0.05)) %>%
      summarize(mean_ratio = mean(ratio),
                mad_ratio = mad(ratio),
                mad3_cutoff = mean_ratio + (3*mad_ratio),
                n_markers = sum(ratio > mad3_cutoff)))


#### Create Ratio plots #### 

load("data/cell_colors.Rdata", verbose = TRUE)

ratio_plots <- map2(marker_stats, names(marker_stats), 
                   ~ggplot(.x, aes(x=ratio, y=std.logFC, color = cellType, shape = rank_ratio <= 5))+
                     geom_point() +
                     facet_wrap(~cellType.target, scales = "free_x") +
                     labs(x = "mean(target logcount)/mean(highest non-target logcount)",
                          title = paste(.y,"cell types"))+ 
                     scale_color_manual(values = cell_colors) +
                     theme_bw()
)
walk2(ratio_plots, names(ratio_plots), ~ggsave(.x, filename = paste0("plots/expr/ratio_vs_stdFC-",.y,".png"), width = 10))


#### Plot expression ####
## rank

exp_plots_ratio <- map2(marker_stats, names(marker_stats), function(x,y){
  ## select right data
  cells <- unique(x$cellType.target)
  names(cells) <- cells
  cellType_col <- ct[[ss(y,"_",2)]]
  temp_sce <- sce[[ss(y,"_")]]
  
  plots <- map(cells,function(cell){
      title = paste(y,cell, "top", top_n ,"markers: Ranked by mean ratios")
      message(title)
      marker_stat_cell <- x %>% filter(cellType.target == cell, rank_ratio <= top_n) %>%
        mutate(Feature = feature_ratio)
      temp_sce_cell <- temp_sce[marker_stat_cell$gene,]
      rownames(temp_sce_cell) <- marker_stat_cell$Feature

      pe <- plotExpression(temp_sce_cell, exprs_values = "logcounts", features = marker_stat_cell$Feature,
                           x=cellType_col, colour_by=cellType_col, point_alpha=0.5, point_size=.2, ncol=5,
                           add_legend=F) +
        scale_color_manual(values = cell_colors)+
        stat_summary(fun = mean, fun.min = mean, fun.max = mean,
                     geom = "crossbar", width = 0.3) +
        geom_text(data = marker_stat_cell, aes(x = -Inf, y = Inf, label = anno),
                  vjust = "inward", hjust = "inward",size = 2.5)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle(label=title)
      ggsave(plot = pe, filename = paste0("plots/expr/",y,"_",cell,"_top",top_n, "_expression_ratio.png"), width = 9, height = 3.25)
      return(pe)
    })
    return(plots)
  })


exp_plots_marker <- map2(marker_stats, names(marker_stats), function(x,y){
  ## select right data
  cells <- unique(x$cellType.target)
  names(cells) <- cells
  cellType_col <- ct[[ss(y,"_",2)]]
  temp_sce <- sce[[ss(y,"_")]]
  
  plots <- map(cells,function(cell){
    title = paste(y,cell, "top", top_n ,"markers: Ranked by stdFC")
    message(title)
    marker_stat_cell <- x %>% filter(cellType.target == cell, rank_marker <= top_n) %>%
      mutate(Feature = feature_marker)
    temp_sce_cell <- temp_sce[marker_stat_cell$gene,]
    rownames(temp_sce_cell) <- marker_stat_cell$Feature
    
    pe <- plotExpression(temp_sce_cell, exprs_values = "logcounts", features = marker_stat_cell$Feature,
                         x=cellType_col, colour_by=cellType_col, point_alpha=0.5, point_size=.2, ncol=5,
                         add_legend=F) +
      scale_color_manual(values = cell_colors)+
      stat_summary(fun = mean, fun.min = mean, fun.max = mean,
                   geom = "crossbar", width = 0.3) +
      geom_text(data = marker_stat_cell, aes(x = -Inf, y = Inf, label = anno),
                vjust = "inward", hjust = "inward",size = 2.5)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ggtitle(label=title)
    ggsave(plot = pe, filename = paste0("plots/expr/",y,"_",cell,"_top",top_n, "expression_marker.png"), width = 9, height = 3.25)
    return(pe)
  })
  return(plots)
})


# sgejobs::job_single('find_markers.R', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript find_markers.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessionsession_info()



