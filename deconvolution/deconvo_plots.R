library(SummarizedExperiment)
library(RColorBrewer)
library(jaffelab)
library(here)
# library(dplyr)
library(reshape2)
# library(ggplot2)
library(patchwork)
library(purrr)
library(tidyverse)
## Load theme colors
source(here("main_colors.R"))

#### Load Data ####
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
pd <- as.data.frame(colData(rse_gene))

## Bind with qSV table
load(here("differential_expression" ,"qSV_mat.Rdata"), verbose = TRUE)
all(rownames(pd)==rownames(qSV_mat))
pd <- cbind(pd, qSV_mat)

## split by region
pd_sacc <- pd[pd$BrainRegion == "sACC",]
pd_amyg <- pd[pd$BrainRegion == "Amygdala",]


## bind with prop data
load("prop_sacc.Rdata", verbose = TRUE)
load("prop_amyg.Rdata", verbose = TRUE)

cells_sacc <- colnames(est_prop_sacc$Est.prop.weighted)
cells_amyg <- colnames(est_prop_amyg$Est.prop.weighted)

## Separate large from small proportions to make plots clear
big_sacc <- map_lgl(est_prop_sacc$Est.prop.weighted, ~mean(.x) > 0.1)
big_amyg <- map_lgl(est_prop_amyg$Est.prop.weighted, ~mean(.x) > 0.1)

prop_sacc <- melt(est_prop_sacc$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop = value)

prop_amyg <- melt(est_prop_amyg$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop = value)

## create long data
#### Compare with top 40 results ####
load("prop_top40_sacc.Rdata", verbose = TRUE)
load("prop_top40_amyg.Rdata", verbose = TRUE)

prop_sacc <- melt(est_prop_top40_sacc$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop_top40 = value) %>%
  right_join(prop_sacc,by = c("sample", "cell_type"))

prop_amyg <- melt(est_prop_top40_amyg$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop_top40 = value)%>%
  right_join(prop_amyg,by = c("sample", "cell_type"))

## Plot
top40_scatter_sacc <- prop40_sacc %>% 
  ggplot(aes(prop, prop_top40, color = cell_type)) +
  geom_point()+ 
  # geom_smooth(method='lm')+
  labs(title = "Prop top 40 vs. all genes, sACC")

ggsave(plot = top40_scatter_sacc , filename = "top40_scatter_sacc.png")

top40_scatter_amyg <- prop40_amyg %>% 
  ggplot(aes(prop, prop_top40, color = cell_type)) +
  geom_point()  + 
  # geom_smooth(method='lm')+
  labs(title = "Prop top 40 vs. all genes, Amygdala")

ggsave(plot = top40_scatter_amyg , filename = "top40_scatter_amyg.png")



#### cell type boxplot ####
prop_dx_sacc <- prop_sacc %>% 
  left_join(pd_sacc %>% select(PrimaryDx) %>% rownames_to_column("sample")) 

prop_dx_amyg <-  prop_amyg %>%
  left_join(pd_amyg %>% select(PrimaryDx) %>% rownames_to_column("sample")) 

## sACC boxplot
bp_dx_sacc_big <- prop_dx_sacc %>% 
  filter(cell_type %in% cells_sacc[big_sacc]) %>%
  ggplot(aes(cell_type, prop, fill = `PrimaryDx`)) +
  geom_boxplot() +
  scale_fill_manual(values = mdd_Dx_colors)+
  labs(title = "Distribution of Cell Type RNA-Propotions",
       subtitle = "sACC samples - MuSiC defaults")+
  theme(legend.position = "None")

bp_dx_sacc_little <- prop_dx_sacc %>% 
  filter(cell_type %in% cells_sacc[!big_sacc]) %>%
  ggplot(aes(cell_type, Prop, fill = `PrimaryDx`)) +
  geom_boxplot() +
  scale_fill_manual(values = mdd_Dx_colors)

ggsave(filename = "cellType_boxplot_sACC.png", 
       plot = bp_dx_sacc_big + bp_dx_sacc_little, 
       width = 10)

## Amyg boxplot
bp_dx_amyg_big <- prop_dx_amyg %>% 
  filter(cell_type %in% cells_amyg[big_amyg]) %>%
  ggplot(aes(cell_type, Prop, fill = `PrimaryDx`)) +
  geom_boxplot() +
  scale_fill_manual(values = mdd_Dx_colors) +
  labs(title = "Distribution of Cell Type RNA-Propotions",
       subtitle = "Amygdala samples - MuSiC defaults") +
  theme(legend.position = "None")

bp_dx_amyg_little <- prop_dx_amyg %>% 
  filter(cell_type %in% cells_amyg[!big_amyg]) %>%
  ggplot(aes(cell_type, Prop, fill = `PrimaryDx`)) +
  geom_boxplot() +
  scale_fill_manual(values = mdd_Dx_colors)

ggsave(filename = "cellType_boxplot_amyg.png", 
       plot = bp_dx_amyg_big + bp_dx_amyg_little, width = 10)

## Top 40 boxplots
## sACC boxplot
bp_dx_sacc_big <- prop_dx_sacc %>% 
  filter(cell_type %in% cells_sacc[big_sacc]) %>%
  ggplot(aes(cell_type, prop_top40, fill = `PrimaryDx`)) +
  geom_boxplot() +
  scale_fill_manual(values = mdd_Dx_colors)+
  labs(title = "Distribution of Cell Type RNA-Propotions",
       subtitle = "sACC samples - top 40 genes")+
  theme(legend.position = "None")

bp_dx_sacc_little <- prop_dx_sacc %>% 
  filter(cell_type %in% cells_sacc[!big_sacc]) %>%
  ggplot(aes(cell_type, prop_top40, fill = `PrimaryDx`)) +
  geom_boxplot() +
  scale_fill_manual(values = mdd_Dx_colors)

ggsave(filename = "cellType_boxplot_top40_sACC.png", 
       plot = bp_dx_sacc_big + bp_dx_sacc_little, 
       width = 10)

## Amyg boxplot
bp_dx_amyg_big <- prop_dx_amyg %>% 
  filter(cell_type %in% cells_amyg[big_amyg]) %>%
  ggplot(aes(cell_type, prop_top40, fill = `PrimaryDx`)) +
  geom_boxplot() +
  scale_fill_manual(values = mdd_Dx_colors) +
  labs(title = "Distribution of Cell Type RNA-Propotions",
       subtitle = "Amygdala samples - top 40 genes") +
  theme(legend.position = "None")

bp_dx_amyg_little <- prop_dx_amyg %>% 
  filter(cell_type %in% cells_amyg[!big_amyg]) %>%
  ggplot(aes(cell_type, prop_top40, fill = `PrimaryDx`)) +
  geom_boxplot() +
  scale_fill_manual(values = mdd_Dx_colors)

ggsave(filename = "cellType_boxplot_top40_amyg.png", 
       plot = bp_dx_amyg_big + bp_dx_amyg_little, width = 10)


#### Prop vs. Age ####
# pd_age <- as.vector(pd_sacc$AgeDeath)
# pd_y <- as.vector(pd_sacc$Oligo)
# png("age_vs_Oligio.png")
# agePlotter(age = pd_age, y = pd_y,
#            mainText = "Avg vs. Prop",
#            mod = model.matrix(~ pd_age),
#            ageBreaks = c(-1, 0, 1, 10, 20, 50, 100))
# dev.off()

age_scatter_sacc <- pd_sacc %>%
  select(all_of(c("RNum", "AgeDeath", cells_sacc))) %>%
  melt(id.vars= c("RNum","AgeDeath")) %>%
  rename(cell_type = variable, prop = value)%>%
  ggplot(aes(AgeDeath, prop, color = cell_type))+
  geom_point()

ggsave(filename = "age_vs_celltype_sACC.png", plot = age_scatter_sacc)

#### check against qSV ####
qsv_names <- colnames(qSV_mat)

prop_vs_qsv_sacc <- pd_sacc %>% select(all_of(c("RNum", cells_sacc))) %>%
  melt(id.vars= c("RNum")) %>%
  rename(cell_type = variable, prop = value) %>%
  left_join(pd_sacc %>% select(all_of(c("RNum", qsv_names))) %>%
              melt(id.vars= c("RNum")) %>%
              rename(qSV_name = variable, qSV = value)) 


scatter_qsv_sacc <- prop_vs_qsv_sacc %>%
  ggplot(aes(qSV, prop)) +
  geom_point(size = .5) +
  facet_grid(cell_type~qSV_name, scales = "free")+
  theme_bw(base_size = 10)

ggsave(filename = "sACC_CellType_qSV.png", plot = scatter_qsv_sacc, width = 26, height = 10)

