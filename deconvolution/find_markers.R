library(SingleCellExperiment)
library(scater)
library(scran)
library(jaffelab)
library(limma)
library(reshape2)
library(tidyverse)
library(patchwork)
library(DeconvoBuddies)
library(here)
library(sessioninfo)

#### Load and organize sce data ####
load(here("deconvolution","data","sce_filtered.Rdata"), verbose = TRUE)
ct <- list(broad = "cellType.Broad", specific = "cellType")

sceXct <- cross2(sce,ct)
names(sceXct) <- map(cross2(names(sce),names(ct)), ~paste0(.x[[1]],"_",.x[[2]]))

#### get marker stats ####
## run get_mean_ratio
mean_ratio <- map(sceXct, ~get_mean_ratio2(sce = .x[[1]], cellType_col = .x[[2]]))
map_int(mean_ratio, nrow)
# sacc_broad    amyg_broad sacc_specific amyg_specific 
# 18581         20330         44783         52444 

## Run findMarkers
markers.t.1vAll <- map(sceXct, ~findMarkers_1vAll(sce = .x[[1]], cellType_col = .x[[2]]))
map_int(markers.t.1vAll, nrow)
# sacc_broad    amyg_broad sacc_specific amyg_specific 
# 106710        106752        177850        213504 

##join stats
marker_stats <- map2(mean_ratio, markers.t.1vAll, ~left_join(.x, .y, by = c("gene", "cellType.target")) %>%
                       mutate(anno = paste0(" ",anno_ratio,"\n",anno_logFC))
                     )

## Save 
save(marker_stats, file = here("deconvolution","data","marker_stats.Rdata"))
# load(here("deconvolution","data", "old","marker_stats.Rdata"), verbose = TRUE)

## save top 5 tables
n_genes <- 5
top_marker_tables <- map2(marker_stats,names(marker_stats), ~filter(.x,rank_ratio<= n_genes) %>% 
                            select(cellType.target,rank_ratio, gene, Symbol, mean.target, cellType.nextHighest = cellType, mean_logcount.nextHighest = mean, std.logFC) %>%
                            arrange(cellType.target) %>%
                            mutate(genecards_url = paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", Symbol),
                                   region = ss(.y,"_"))
)

walk2(top_marker_tables,names(marker_stats), ~write_csv(.x, paste0("data/top",n_genes,"_markers_",.y,".csv")))  

marker_genes <- map(top_marker_tables, ~pull(.x, "gene"))
save(marker_genes, file = here("deconvolution","data","marker_genes.Rdata"))

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

walk2(sceXct, names(sceXct), function(sxc, name_sxc){
  message(name_sxc)
  pdf(here("deconvolution","plots","expr",paste0("expr_top10_mean_ratio-",name_sxc,".pdf")))
  walk(levels(marker_stats[[name_sxc]]$cellType.target),
       
       ~print(plot_marker_express(sxc[[1]], 
                                  marker_stats[[name_sxc]], 
                                  cell_type = .x, 
                                  n_genes = 10, 
                                  rank_col = "rank_ratio", 
                                  anno_col = "anno", 
                                  cellType_col = sxc[[2]])+
                scale_color_manual(values = cell_colors)
       )
  )
  dev.off()
  
})


walk2(sceXct, names(sceXct), function(sxc, name_sxc){
  message(name_sxc)
  pdf(here("deconvolution","plots","expr",paste0("expr_top10_1vAll-",name_sxc,".pdf")))
  walk(levels(marker_stats[[name_sxc]]$cellType.target),
       
       ~print(plot_marker_express(sxc[[1]], 
                                  marker_stats[[name_sxc]], 
                                  cell_type = .x, 
                                  n_genes = 10, 
                                  rank_col = "rank_marker", 
                                  anno_col = "anno", 
                                  cellType_col = sxc[[2]])+
                scale_color_manual(values = cell_colors)
       )
  )
  dev.off()
  
})


# sgejobs::job_single('find_markers', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript find_markers.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



