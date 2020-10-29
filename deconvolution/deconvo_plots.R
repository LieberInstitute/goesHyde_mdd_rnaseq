library(RColorBrewer)
library(jaffelab)
library(here)
library(dplyr)
library(reshape2)

#### plot ####
# load("prop_sacc.Rdata", verbose = TRUE)
pd_sacc <- as.data.frame(colData(rse_gene_sACC))
pd_sacc <- cbind(pd_sacc, est_prop_sacc$Est.prop.weighted)

bp_dx_sacc <- pd_sacc %>% select(PrimaryDx, colnames(est_prop_sacc$Est.prop.weighted)) %>%
  melt(id.vars = "PrimaryDx") %>%
  rename(`Cell Type` = variable, Prop = value) %>%
  ggplot(aes(`Cell Type`, Prop, fill = `PrimaryDx`)) +
  geom_boxplot() +
  labs(title = "Distribution of Cell Type Propotions",
       subtitle = "sACC samples - MuSiC defaults")

ggsave(filename = "cellType_boxplot_sACC.png", plot = bp_dx_sacc, width = 10)

# load("prop_sacc.Rdata", verbose = TRUE)
pd_sacc <- as.data.frame(colData(rse_gene_sACC))
pd_sacc <- cbind(pd_sacc, est_prop_sacc$Est.prop.weighted)

bp_dx_amyg <- pd_sacc %>% select(PrimaryDx, names(cell_bias$M.S)) %>%
  melt(id.vars = "PrimaryDx") %>%
  rename(`Cell Type` = variable, Prop = value) %>%
  ggplot(aes(`Cell Type`, Prop, fill = `PrimaryDx`)) +
  geom_boxplot() +
  labs(title = "Distribution of Cell Type Propotions",
       subtitle = "sACC samples - MuSiC defaults")

ggsave(filename = "sACC_CellType_boxplot.png", plot = bp_sacc_dx, width = 10)

#### check against qSV 1####
load(here("differential_expression" ,"qSV_mat.Rdata"), verbose = TRUE)
all(colnames(rse_gene)==rownames(qSV_mat))

pd_sacc <- cbind(pd_sacc, qSV_mat[rownames(pd_sacc),])

scatter_qsv1 <- pd_sacc %>% select(PrimaryDx, colnames(est_prop_sacc_noT$Est.prop.weighted), PC1) %>%
  melt(id.vars = c("PrimaryDx","PC1")) %>% 
  rename(`Cell Type` = variable, Prop = value) %>%
  ggplot(aes(PC1, Prop, color = `Cell Type`)) +
  geom_point() +
  facet_wrap(~`Cell Type`, scales = "free_y")+
  labs(title = "qSV1 vs. Cell Type Prop",
       subtitle = "sACC samples - MuSiC defaults - no T cell")

ggsave(filename = "sACC_CellType_qSV1.png", plot = scatter_qsv1, width = 14)
