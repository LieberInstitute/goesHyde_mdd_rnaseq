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
source("big_little_boxplot.R")
## extract cell_types establish color pallet
load(here("deconvolution","data","cell_colors.Rdata"), verbose = TRUE)
load(here("deconvolution","data","prop_long.Rdata"), verbose = TRUE)

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
age_scatter_broad_sacc <- prop_age_broad_sacc %>%
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


#### check against qSV ####
qsv_names <- colnames(qSV_mat)

qsv_long_sacc <- pd_sacc %>% 
  rownames_to_column("sample") %>%
  select(all_of(c("sample", qsv_names))) %>%
  melt(id.vars= c("sample")) %>%
  rename(qSV_name = variable, qSV = value) 

qsv_long_amyg <- pd_amyg %>% 
  rownames_to_column("sample") %>%
  select(all_of(c("sample", qsv_names))) %>%
  melt(id.vars= c("sample")) %>%
  rename(qSV_name = variable, qSV = value) 

prop_qsv_sacc <- prop_sacc %>%
  left_join(qsv_long_sacc, by = "sample")

prop_qsv_amyg <- prop_amyg %>%
  left_join(qsv_long_amyg, by = "sample")

#### Calculate p-values ####
##sACC
cross_sacc <- cross2(cells_sacc, qsv_names)

cross_df_sacc <- as.data.frame(do.call("rbind", cross_sacc))
colnames(cross_df_sacc) <- c("cell_type","qSV_name")
cross_df_sacc$cell_type <- factor(cross_df_sacc$cell_type, levels = unique(cross_df_sacc$cell_type))
cross_df_sacc$qSV_name <- factor(cross_df_sacc$qSV_name, levels = unique(cross_df_sacc$qSV_name))

## All genes
qsv_lm_sacc <- bind_rows(map(cross_sacc, function(x){
  cell = x[[1]]
  qsv = x[[2]]
  prop_qsv <- prop_qsv_sacc %>%
    filter(cell_type == cell & qSV_name == qsv) 
  
  tidy(lm(prop~qSV-1, data = prop_qsv))
}))

qsv_lm_sacc_df <- cbind(cross_df_sacc,qsv_lm_sacc) 

qsv_lm_sacc_df$p.bonf <- p.adjust(qsv_lm_sacc_df$p.value, 'bonf')
qsv_lm_sacc_df$p.bonf.sig <- qsv_lm_sacc_df$p.bonf < 0.05
qsv_lm_sacc_df$p.fdr<- p.adjust(qsv_lm_sacc_df$p.value, 'fdr')

sum(qsv_lm_sacc_df$p.bonf < 0.05)
# [1] 74
head(qsv_lm_sacc_df)

## top40 
qsv_lm_top40_sacc <- bind_rows(map(cross_sacc, function(x){
  cell = x[[1]]
  qsv = x[[2]]
  prop_qsv <- prop_qsv_sacc %>%
    filter(cell_type == cell & qSV_name == qsv) 
  
  tidy(lm(prop_top40~qSV-1, data = prop_qsv))
}))

qsv_lm_top40_sacc_df <- cbind(cross_df_sacc,qsv_lm_top40_sacc) 

qsv_lm_top40_sacc_df$p.bonf <- p.adjust(qsv_lm_top40_sacc_df$p.value, 'bonf')
qsv_lm_top40_sacc_df$p.bonf.sig <- qsv_lm_top40_sacc_df$p.bonf < 0.05
qsv_lm_top40_sacc_df$p.fdr<- p.adjust(qsv_lm_top40_sacc_df$p.value, 'fdr')

sum(qsv_lm_top40_sacc_df$p.bonf < 0.05)
# [1] 83
head(qsv_lm_top40_sacc_df)

## Amyg
cross_amyg <- cross2(cells_amyg, qsv_names)

cross_df_amyg <- as.data.frame(do.call("rbind", cross_amyg))
colnames(cross_df_amyg) <- c("cell_type","qSV_name")
cross_df_amyg$cell_type <- factor(cross_df_amyg$cell_type, levels = unique(cross_df_amyg$cell_type))
cross_df_amyg$qSV_name <- factor(cross_df_amyg$qSV_name, levels = unique(cross_df_amyg$qSV_name))

## All genes
qsv_lm_amyg <- bind_rows(map(cross_amyg, function(x){
  cell = x[[1]]
  qsv = x[[2]]
  prop_qsv <- prop_qsv_amyg %>%
    filter(cell_type == cell & qSV_name == qsv) 
  
  tidy(lm(prop~qSV-1, data = prop_qsv))
}))

qsv_lm_amyg_df <- cbind(cross_df_amyg, qsv_lm_amyg) 

qsv_lm_amyg_df$p.bonf <- p.adjust(qsv_lm_amyg_df$p.value, 'bonf')
qsv_lm_amyg_df$p.bonf.sig <- qsv_lm_amyg_df$p.bonf < 0.05
qsv_lm_amyg_df$p.fdr<- p.adjust(qsv_lm_amyg_df$p.value, 'fdr')

sum(is.na(qsv_lm_amyg_df$p.value))
# [1] 26
## Replace NA with FALSE 
qsv_lm_amyg_df$p.bonf.sig[is.na(qsv_lm_amyg_df$p.bonf.sig)] <- FALSE
sum(qsv_lm_amyg_df$p.bonf.sig)
# [1] 69

## top40
qsv_lm_top40_amyg <- bind_rows(map(cross_amyg, function(x){
  cell = x[[1]]
  qsv = x[[2]]
  prop_qsv <- prop_qsv_amyg %>%
    filter(cell_type == cell & qSV_name == qsv) 
  
  tidy(lm(prop_top40~qSV-1, data = prop_qsv))
}))

qsv_lm_top40_amyg_df <- cbind(cross_df_amyg, qsv_lm_top40_amyg) 

qsv_lm_top40_amyg_df$p.bonf <- p.adjust(qsv_lm_top40_amyg_df$p.value, 'bonf')
qsv_lm_top40_amyg_df$p.bonf.sig <- qsv_lm_top40_amyg_df$p.bonf < 0.05
qsv_lm_top40_amyg_df$p.fdr<- p.adjust(qsv_lm_top40_amyg_df$p.value, 'fdr')

sum(qsv_lm_top40_amyg_df$p.bonf.sig)
# [1] 85

#### Create scatter plots ####
scatter_qsv_sacc <- prop_qsv_sacc %>% 
  left_join(qsv_lm_sacc_df %>% select(cell_type, qSV_name, p.bonf.sig)) %>% 
  ggplot(aes(qSV, prop, color= p.bonf.sig)) +
  geom_point(size = .4, alpha = .25) +
  facet_grid(cell_type~qSV_name, scales = "free")+
  theme_bw(base_size = 10)+
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'red'))

ggsave(filename = "plots/cellType_qSV_sACC.png", plot = scatter_qsv_sacc, width = 26, height = 10)

scatter_qsv_top40_sacc <- prop_qsv_sacc %>% 
  left_join(qsv_lm_top40_sacc_df %>% select(cell_type, qSV_name, p.bonf.sig)) %>% 
  ggplot(aes(qSV, prop_top40, color= p.bonf.sig)) +
  geom_point(size = .4, alpha = .25) +
  facet_grid(cell_type~qSV_name, scales = "free")+
  theme_bw(base_size = 10)+
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'red'))

ggsave(filename = "plots/cellType_qSV_top40_sACC.png", plot = scatter_qsv_top40_sacc, width = 26, height = 10)

scatter_qsv_amyg <- prop_qsv_amyg %>% 
  left_join(qsv_lm_amyg_df %>% select(cell_type, qSV_name, p.bonf.sig)) %>% 
  ggplot(aes(qSV, prop, color= p.bonf.sig)) +
  geom_point(size = .5, alpha = .25) +
  facet_grid(cell_type~qSV_name, scales = "free")+
  theme_bw(base_size = 10)+
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'red'))

ggsave(filename = "plots/cellType_qSV_amyg.png", plot = scatter_qsv_amyg, width = 26, height = 10)

scatter_qsv_top40_amyg <- prop_qsv_amyg %>% 
  left_join(qsv_lm_top40_amyg_df %>% select(cell_type, qSV_name, p.bonf.sig)) %>% 
  ggplot(aes(qSV, prop_top40, color= p.bonf.sig)) +
  geom_point(size = .5, alpha = .25) +
  facet_grid(cell_type~qSV_name, scales = "free")+
  theme_bw(base_size = 10)+
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'red'))

ggsave(filename = "plots/cellType_qSV_top40_amyg.png", plot = scatter_qsv_top40_amyg, width = 26, height = 10)

# sgejobs::job_single('deconvo_plots', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript deconvo_plots.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessionsession_info()
