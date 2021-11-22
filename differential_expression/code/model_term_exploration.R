library(SummarizedExperiment)
library(tidyverse)
library(here)
library(sessioninfo)

load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
load(here("data", "MDD_colors.Rdata"), verbose = TRUE)
load(here("differential_expression","data","qSV_mat.Rdata"), verbose = TRUE)


pd <- as.data.frame(colData(rse_gene))
pd_long <- pd %>% select(PrimaryDx, Experiment, BrainRegion, 
                         #AgeDeath,  snpPC1 ,snpPC2, snpPC3, snpPC4, snpPC5,
                           mitoRate, rRNA_rate,totalAssignedGene, RIN, ERCCsumLogErr) %>%
  mutate(ERCCsumLogErr_abs = abs(ERCCsumLogErr)) %>%
  pivot_longer(!c(PrimaryDx, Experiment, BrainRegion), names_to = "Metric")


qc_boxplots <- ggplot(pd_long, aes(x = BrainRegion, y = value, fill = PrimaryDx, color = Experiment)) +
  geom_boxplot() +
  facet_wrap(~Metric, scales = "free_y", nrow = 2) +
  scale_color_manual(values = c(psychENCODE_MDD = "black", psychENCODE_BP = "skyblue")) +
  scale_fill_manual(values = mdd_Dx_colors) +
  theme_bw() +
  theme(text = element_text(size=15)) 

ggsave(qc_boxplots, file = here("differential_expression","plots","model_terms_boxplots.png"), width = 15)

#### qSV plots ####
qSV_long <- pd %>% select(PrimaryDx, Experiment, BrainRegion) %>%
  cbind(qSV_mat) %>%
  pivot_longer(!c(PrimaryDx, Experiment, BrainRegion), names_to = "qSV")

qSV_long$qSV <- factor(qSV_long$qSV, levels = paste0("PC", 1:26))
levels(qSV_long$qSV)

qSV_boxplots <- ggplot(qSV_long, aes(x = BrainRegion, y = value, fill = PrimaryDx, color = Experiment)) +
  geom_boxplot() +
  facet_wrap(~qSV, scales = "free_y", nrow = 3) +
  scale_color_manual(values = c(psychENCODE_MDD = "black", psychENCODE_BP = "skyblue")) +
  scale_fill_manual(values = mdd_Dx_colors) +
  theme_bw() +
  theme(text = element_text(size=15), legend.position="top") 

ggsave(qSV_boxplots, file = here("differential_expression","plots","model_qSV_boxplots.png"), width = 20, height = 15)

