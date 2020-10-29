library(SummarizedExperiment)
library(RColorBrewer)
library(jaffelab)
library(here)
library(dplyr)
library(reshape2)
library(ggplot2)

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

source(here("main_colors.R"))

## bind with prop data
load("prop_sacc.Rdata", verbose = TRUE)
load("prop_amyg.Rdata", verbose = TRUE)

pd_sacc <- cbind(pd_sacc, est_prop_sacc$Est.prop.weighted)
pd_amyg <- cbind(pd_amyg, est_prop_amyg$Est.prop.weighted)

cells_sacc <- colnames(est_prop_sacc$Est.prop.weighted)
cells_amyg <- colnames(est_prop_amyg$Est.prop.weighted)

#### cell type boxplot ####

bp_dx_sacc <- pd_sacc %>% select(PrimaryDx, cells_sacc) %>%
  melt(id.vars = "PrimaryDx") %>%
  rename(`Cell Type` = variable, Prop = value) %>%
  ggplot(aes(`Cell Type`, Prop, fill = `PrimaryDx`)) +
  geom_boxplot() +
  scale_fill_manual(values = mdd_Dx_colors)+
  labs(title = "Distribution of Cell Type Propotions",
       subtitle = "sACC samples - MuSiC defaults")

ggsave(filename = "cellType_boxplot_sACC.png", plot = bp_dx_sacc, width = 10)


bp_dx_amyg <- pd_amyg %>% select(all_of(c("PrimaryDx", cells_amyg))) %>%
  melt(id.vars = "PrimaryDx") %>%
  rename(`Cell Type` = variable, Prop = value) %>%
  ggplot(aes(`Cell Type`, Prop, fill = `PrimaryDx`)) +
  geom_boxplot() +
  scale_fill_manual(values = mdd_Dx_colors) +
  labs(title = "Distribution of Cell Type Propotions",
       subtitle = "Amygdala samples - MuSiC defaults")

ggsave(filename = "cellType_boxplot_amyg.png", plot = bp_dx_amyg, width = 10)

#### check against qSV 1####

scatter_qsv1 <- pd_sacc %>% select(PrimaryDx, colnames(est_prop_sacc_noT$Est.prop.weighted), PC1) %>%
  melt(id.vars = c("PrimaryDx","PC1")) %>% 
  rename(`Cell Type` = variable, Prop = value) %>%
  ggplot(aes(PC1, Prop, color = `Cell Type`)) +
  geom_point() +
  facet_wrap(~`Cell Type`, scales = "free_y")+
  labs(title = "qSV1 vs. Cell Type Prop",
       subtitle = "sACC samples - MuSiC defaults - no T cell")

ggsave(filename = "sACC_CellType_qSV1.png", plot = scatter_qsv1, width = 14)
