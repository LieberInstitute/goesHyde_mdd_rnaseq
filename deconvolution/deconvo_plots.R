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

source(here("main_colors.R"))
source(here("deconvolution","big_little_boxplot.R"))
## extract cell_types establish color pallet
load(here("deconvolution","data","cell_colors.Rdata"), verbose = TRUE)
# load(here("deconvolution","data","prop_long.Rdata"), verbose = TRUE)
load(here("deconvolution","data","est_prop_top5_long.Rdata"), verbose = TRUE)


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

cells_broad <- c("Inhib", "Excit","Astro","Micro","Oligo","OPC")
## Plot
top40_scatter_sacc <- prop_sacc %>% 
  ggplot(aes(prop, prop_top40, color = cell_type)) +
  geom_point(size = 0.5)+ 
  labs(title = "Prop top-40 vs. all genes",
       subtitle = "sACC")+
  geom_abline()+
  facet_wrap(~cell_type, scales = "free")+
  scale_color_manual(values = cell_colors)+
  theme(legend.position = "None")

top40_scatter_amyg <- prop_amyg %>% 
  ggplot(aes(prop, prop_top40, color = cell_type)) +
  geom_point(size = 0.5)  + 
  labs(subtitle = "Amygdala") +
  geom_abline()+
  facet_wrap(~cell_type, scales = "free")+
  scale_color_manual(values = cell_colors)+
  theme(legend.position = "None")

ggsave(plot = top40_scatter_sacc + top40_scatter_amyg , filename = "plots/top40_scatter.png", width = 16)

## Plot broad
top40_scatter_broad_sacc <- prop_broad_sacc %>% 
  ggplot(aes(prop, prop_top40, color = cell_type)) +
  geom_point(size = 0.5)+ 
  geom_abline()+
  scale_color_manual(values = cell_colors)+
  labs(title = "Prop top-40 vs. all genes",
       subtitle = "sACC")+
  facet_wrap(~cell_type, scales = "free")+
  theme(legend.position = "None")

top40_scatter_broad_amyg <- prop_broad_amyg %>% 
  ggplot(aes(prop, prop_top40, color = cell_type)) +
  geom_point(size = 0.5)  +  
  geom_abline()+
  scale_color_manual(values = cell_colors)+
  labs(subtitle = "Amygdala") +
  facet_wrap(~cell_type, scales = "free")+
  theme(legend.position = "None")

ggsave(plot = top40_scatter_broad_sacc + top40_scatter_broad_amyg , filename = "plots/top40_scatter_broad.png", width = 16)

#### cell type boxplot ####

## specific cell types
both_regions_bl_boxplot(prop_sacc, prop_amyg, yvar = "prop", 
                        title = "Cell Type RNA-prop - all common genes",
                        filename = "plots/cellType_boxplots.png")

both_regions_bl_boxplot(prop_sacc, prop_amyg, yvar = "prop_top40", 
                        title = "Cell Type RNA-prop - top 40 genes",
                        filename = "plots/cellType_boxplots_top40.png")

## broad cell types
both_regions_bl_boxplot(prop_broad_sacc, prop_broad_amyg, yvar = "prop", 
                        title = "Cell Type RNA-prop -  all common genes", 
                        filename = "plots/cellType_boxplots_broad.png")

both_regions_bl_boxplot(prop_broad_sacc, prop_broad_amyg, yvar = "prop_top40", 
                        title = "Cell Type RNA-prop - top 40 genes", 
                        filename = "plots/cellType_boxplots_broad_top40.png")

## Age scatter
age_scatter_broad_sacc <- prop_broad_sacc %>%
  ggplot(aes(AgeDeath, prop, color = cell_type))+
  geom_point() +
  facet_wrap(~cell_type, scales = "free_y", nrow = 6)

ggsave(filename = "plots/age_vs_celltype_broad_sACC.png", plot = age_scatter_broad_sacc)

#### Compare summed specific cell types with broad ####
prop_sum_sacc <- prop_sacc %>% 
  filter(!cell_type %in% cells_broad) %>%
  mutate(cell_type = gsub(".[0-9]", "", cell_type)) %>%
  group_by(sample, cell_type) %>%
  summarise(sum_prop = sum(prop),
            sum_prop_top40 = sum(prop_top40)) %>%
  left_join(prop_broad_sacc, by = c("sample","cell_type"))

prop_sum_amyg <- prop_amyg %>% 
  filter(!cell_type %in% cells_broad) %>%
  mutate(cell_type = gsub(".[0-9]", "", cell_type)) %>%
  group_by(sample, cell_type) %>%
  summarise(sum_prop = sum(prop),
            sum_prop_top40 = sum(prop_top40)) %>%
  left_join(prop_broad_amyg, by = c("sample","cell_type"))

## plot
sum_scatter_sacc <- prop_sum_sacc %>%
  ggplot(aes(prop, sum_prop, color = cell_type)) +
  geom_point()  +  
  geom_abline()+
  scale_color_manual(values = cell_colors)+
  facet_wrap(~cell_type, scales = "free")+
  theme(legend.position = "None") +
  labs(title = "Specific prop summed vs. Broad prop",
       subtitle = "all genes - sACC")

sum_scatter_amyg <- prop_sum_amyg %>%
  ggplot(aes(prop, sum_prop, color = cell_type)) +
  geom_point()  +  
  geom_abline()+
  scale_color_manual(values = cell_colors)+
  facet_wrap(~cell_type, scales = "free")+
  theme(legend.position = "None") +
  labs(subtitle = "Amygdala")

ggsave(filename = "plots/sum_scatter.png",
       plot = sum_scatter_sacc + sum_scatter_amyg, width = 14)

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
