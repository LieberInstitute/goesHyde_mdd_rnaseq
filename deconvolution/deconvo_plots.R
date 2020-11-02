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


## bind with prop data
load("prop_sacc.Rdata", verbose = TRUE)
load("prop_amyg.Rdata", verbose = TRUE)

cells_sacc <- colnames(est_prop_sacc$Est.prop.weighted)
cells_amyg <- colnames(est_prop_amyg$Est.prop.weighted)

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
## Add dx to long data
prop_dx_sacc <- prop_sacc %>% 
  left_join(pd_sacc %>% select(PrimaryDx) %>% rownames_to_column("sample")) 

prop_dx_amyg <-  prop_amyg %>%
  left_join(pd_amyg %>% select(PrimaryDx) %>% rownames_to_column("sample")) 

## sACC boxplot
big_little_boxplot(prop_dx_sacc, xvar = "cell_type", yvar = "prop",
                   fillvar =  "PrimaryDx",
                   title = "Distribution of Cell Type RNA-Propotions",
                   subtitle = "sACC samples - MuSiC defaults", 
                   fn = "cellType_boxplot_sACC.png")

## Amyg boxplot
big_little_boxplot(prop_dx_amyg, xvar = "cell_type", yvar = "prop",
                   fillvar =  "PrimaryDx",
                   title = "Distribution of Cell Type RNA-Propotions",
                   subtitle = "Amygdala samples - MuSiC defaults", 
                   fn = "cellType_boxplot_amyg.png")

## Top 40 boxplots
## sACC boxplot
big_little_boxplot(prop_dx_sacc, xvar = "cell_type", yvar = "prop_top40",
                   fillvar =  "PrimaryDx",
                   title = "Distribution of Cell Type RNA-Propotions",
                   subtitle = "sACC samples - top 40 genes", 
                   fn = "cellType_boxplot_top40_sACC.png")

## Amyg boxplot
big_little_boxplot(prop_dx_amyg, xvar = "cell_type", yvar = "prop_top40",
                   fillvar =  "PrimaryDx",
                   title = "Distribution of Cell Type RNA-Propotions",
                   subtitle = "Amygdala samples - top 40 genes", 
                   fn = "cellType_boxplot_top40_amyg.png")


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


