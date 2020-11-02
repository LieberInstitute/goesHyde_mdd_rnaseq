library(SummarizedExperiment)
library(RColorBrewer)
library(jaffelab)
library(here)
library(reshape2)
library(patchwork)
library(purrr)
library(tidyverse)

source(here("main_colors.R"))
source("big_little_boxplot.R")
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


## load and melt prop data
load("prop_sacc.Rdata", verbose = TRUE)
load("prop_amyg.Rdata", verbose = TRUE)
cells_sacc <- colnames(est_prop_sacc$Est.prop.weighted)
cells_amyg <- colnames(est_prop_amyg$Est.prop.weighted)

prop_sacc <- melt(est_prop_sacc$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop = value)

prop_amyg <- melt(est_prop_amyg$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop = value)

## broad data
load("prop_broad_sacc.Rdata", verbose = TRUE)
load("prop_broad_amyg.Rdata", verbose = TRUE)

prop_broad_sacc <- melt(est_prop_broad_sacc$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop = value) 

prop_broad_amyg <- melt(est_prop_broad_amyg$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop = value) 

#### Compare with top 40 results ####
load("prop_top40_sacc.Rdata", verbose = TRUE)
load("prop_top40_amyg.Rdata", verbose = TRUE)

prop_sacc <- melt(est_prop_top40_sacc$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop_top40 = value) %>%
  right_join(prop_sacc,by = c("sample", "cell_type"))

prop_amyg <- melt(est_prop_top40_amyg$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop_top40 = value)%>%
  right_join(prop_amyg,by = c("sample", "cell_type"))

load("prop_top40_broad_sacc.Rdata", verbose = TRUE)
load("prop_top40_broad_amyg.Rdata", verbose = TRUE)

prop_broad_sacc <- melt(est_prop_top40_broad_sacc$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop_top40 = value) %>%
  right_join(prop_broad_sacc,by = c("sample", "cell_type"))

prop_broad_amyg <- melt(est_prop_top40_broad_amyg$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop_top40 = value)%>%
  right_join(prop_broad_amyg,by = c("sample", "cell_type"))


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
## Add dx to long data
prop_dx_sacc <- prop_sacc %>% 
  left_join(pd_sacc %>% select(PrimaryDx) %>% rownames_to_column("sample")) 

prop_broad_dx_sacc <- prop_broad_sacc %>% 
  left_join(pd_sacc %>% select(PrimaryDx) %>% rownames_to_column("sample")) 

prop_broad_dx_amyg <-  prop_broad_amyg %>%
  left_join(pd_amyg %>% select(PrimaryDx) %>% rownames_to_column("sample")) 

## specific cell types
both_regions_bl_boxplot(prop_dx_sacc, prop_dx_amyg, yvar = "prop", 
                        title = "Cell Type RNA-prop - all common genes",
                        filename = "cellType_boxplots.png")

both_regions_bl_boxplot(prop_dx_sacc, prop_dx_amyg, yvar = "prop_top40", 
                        title = "Cell Type RNA-prop - all common genes",
                        filename = "cellType_boxplots_top40.png")

## broad cell types
both_regions_bl_boxplot(prop_broad_dx_sacc, prop_broad_dx_amyg, yvar = "prop", 
                        title = "Cell Type RNA-prop - top 40 genes", 
                        filename = "cellType_boxplots_broad.png")

both_regions_bl_boxplot(prop_broad_dx_sacc, prop_broad_dx_amyg, yvar = "prop_top40", 
                        title = "Cell Type RNA-prop - top 40 genes", 
                        filename = "cellType_boxplots_broad_top40.png")



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


