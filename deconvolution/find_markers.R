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
load(here("deconvolution","data", "cell_colors.Rdata"), verbose = TRUE)

#### get marker stats ####
## run get_mean_ratio
mean_ratio <- get_mean_ratio2(sce = sce_pan, cellType_col = "cellType.Broad")
nrow(mean_ratio)
# [1] 17179

## Run findMarkers
markers.t.1vAll <- findMarkers_1vAll(sce = sce_pan, cellType_col = "cellType.Broad")
nrow(markers.t.1vAll)

##join stats
marker_stats <- left_join(mean_ratio, markers.t.1vAll, by = c("gene", "cellType.target")) %>%
  mutate(anno = paste0(" ",anno_ratio,"\n",anno_logFC))

## Save 
save(marker_stats, file = here("deconvolution","data","marker_stats.Rdata"))

## save marker gene tables
n_genes <- 25
top_marker_table <- marker_stats %>% 
  filter(rank_ratio <= n_genes) %>% 
  select(cellType.target,rank_ratio, gene, Symbol, mean.target, cellType.nextHighest = cellType, mean_logcount.nextHighest = mean, std.logFC) %>%
  arrange(cellType.target) %>%
  mutate(genecards_url = paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", Symbol))

top_marker_table %>% write_csv(here("deconvolution","data",paste0("top",n_genes,"_markers.csv")))

marker_genes <- top_marker_table %>% pull("gene")
save(marker_genes, file = here("deconvolution","data","marker_genes.Rdata"))

#### Create Ratio plots #### 
ratio_plot <- marker_stats %>% mutate(Marker = rank_ratio <= n_genes) %>%
  ggplot(aes(x=ratio, y=std.logFC, color = Marker))+
  geom_point() +
  facet_wrap(~cellType.target, scales = "free_x") +
  labs(x = "mean(target logcount)/mean(highest non-target logcount)")+ 
  # scale_color_manual(values = cell_colors) +
  theme_bw()

ggsave(ratio_plot, filename = paste0("plots/expr/ratio_vs_stdFC.png"), width = 10)

#### Plot expression ####
# pdf(here("deconvolution","plots","expr","expr_mean_ratio.pdf"))
walk(levels(marker_stats$cellType.target)[1],
     # ~message(.x)
     ~print(DeconvoBuddies::plot_marker_express(sce = sce_pan,
                                stats = marker_stats,
                                cell_type = .x,
                                n_genes = 25,
                                rank_col = "rank_ratio",
                                anno_col = "anno",
                                cellType_col = "cellType.Broad"+
              scale_color_manual(values = cell_colors)
     )
)
)
dev.off()

pdf(here("deconvolution","plots","expr","expr_mean_ratio.pdf"))
walk(levels(marker_stats$cellType.target[1]),
     
     ~print(plot_marker_express(sce_pan, 
                                marker_stats, 
                                cell_type = .x, 
                                n_genes = 25, 
                                rank_col = "rank_marker", 
                                anno_col = "anno", 
                                cellType_col = "cellType.Broad"+
                                  scale_color_manual(values = cell_colors)
     ))
)
dev.off()

# sgejobs::job_single('find_markers', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript find_markers.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



