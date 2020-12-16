library(SummarizedExperiment)
library(RColorBrewer)
library(jaffelab)
library(here)
library(reshape2)
library(patchwork)
library(purrr)
library(tidyverse)
library(broom)
library(sessioninfo)
library(pheatmap)

## Load colors and plotting functions
source(here("main_colors.R"))
source(here("deconvolution","big_little_boxplot.R"))
load(here("deconvolution","data","cell_colors.Rdata"), verbose = TRUE)

## Load MuSiC results & res data
load(here("deconvolution","data","est_prop_top5_long.Rdata"), verbose = TRUE)
load(here("deconvolution","data","est_prop_top5.Rdata"),verbose = TRUE)
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)

#### Boxplots ####
walk2(est_prop_long, names(est_prop_long), function(x,y){
  blb <- big_little_boxplot(data = x,
                            xvar = "cell_type", 
                            yvar = "prop",
                            fillvar =  "PrimaryDx",
                            colorvar = "ignore",
                            title = "MuSiC Proptions: Top5 markers",
                            subtitle = y)
  ggsave(plot = blb, filename = paste0("plots/cellType_boxplots_",y,".png"), width = 15)
} )

walk2(est_prop_long, names(est_prop_long), function(x,y){
  blb <- big_little_boxplot(data = x,
                            xvar = "cell_type", 
                            yvar = "prop",
                            fillvar =  "PrimaryDx",
                            colorvar = "Sex",
                            title = "MuSiC Proptions: Top5 markers",
                            subtitle = y)
  ggsave(plot = blb, filename = paste0("plots/cellType_sex_boxplots_",y,".png"), width = 15)
} )

#### Experimental plots ####
bar <-est_prop_long[[1]]%>%
  left_join(est_prop_long[[1]] %>%
              filter(cell_type == 'Oligo') %>%
              select(sample, prop_1 = prop),
            by = 'sample') %>% 
  ggplot(aes(x = reorder(sample, desc(prop_1)), y = prop, fill = cell_type)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_colors)+
  coord_flip()

ggsave(bar, filename = "plots/bar_test.png", height = 10)

## try pheatmap
my_anno_colors <- list(cellType = cell_colors,
                       PrimaryDx = mdd_Dx_colors,
                       Experiment = mdd_dataset_colors)

sample_anno <- as.data.frame(colData(rse_gene)) %>%
  select(Experiment, Sex, PrimaryDx, AgeDeath)

walk2(est_prop, names(est_prop), function(est, n){
  epw <- est$Est.prop.weighted
  png(paste0("plots/MuSiC-heatmap_",n,".png"), height = 800, width = 580)
  pheatmap(epw,
           show_rownames = FALSE,
           annotation_row = sample_anno,
           annotation_colors = my_anno_colors,
           main = paste(n,"cell type proptions"))
  dev.off()
})




## plot top 40 
sum_scatter_sacc <- prop_sum_sacc %>%
  ggplot(aes(prop_top40, sum_prop_top40, color = cell_type)) +
  geom_point()  +  
  geom_abline()+
  scale_color_manual(values = cell_colors)+
  facet_wrap(~cell_type, scales = "free")+
  theme(legend.position = "None") +
  labs(title = "Specific prop summed vs. Broad prop",
       subtitle = "top 40 - sACC")

sum_scatter_amyg <- prop_sum_amyg %>%
  ggplot(aes(prop_top40, sum_prop_top40, color = cell_type)) +
  geom_point()  +  
  geom_abline()+
  scale_color_manual(values = cell_colors)+
  facet_wrap(~cell_type, scales = "free")+
  theme(legend.position = "None") +
  labs(subtitle = "Amygdala")

ggsave(filename = "plots/sum_scatter_top40.png",
       plot = sum_scatter_sacc + sum_scatter_amyg, width = 14)


# sgejobs::job_single('deconvo_plots', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript deconvo_plots.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessionsession_info()
