library(SummarizedExperiment)
library(tidyverse)
library(broom)
library(here)
library(sessioninfo)

load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
load(here("data", "MDD_colors.Rdata"), verbose = TRUE)
load(here("differential_expression","data","qSV_mat.Rdata"), verbose = TRUE)


pd <- as.data.frame(colData(rse_gene))

#### QC metrics ####
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

#### qSVs vs. QC metrics ####
pd_qsv <- pd_long %>% 
  rename(metric_value = value) %>%
  left_join(qSV_long %>% rename(qsv_value = value))

pd_qsv_scatter <- pd_qsv %>%
  ggplot(aes(qSV_value, metric_value, color = Experiment)) +
  geom_point() +
  facet_grid(Metric~qSV)

ggsave(pd_qsv_scatter, file = here("differential_expression","plots","model_qSV_metric_scatter.png"), width = 20, height = 15)

pd_qSV_fit <- pd_qsv %>% group_by(Metric, qSV, BrainRegion) %>%
  do(fitQSV = tidy(lm(metric_value ~ qSV_value, data = .))) %>%
  unnest(fitQSV) %>%
  filter(term == "qSV_value") %>%
  mutate(p.bonf = p.adjust(p.value, "bonf"),
         p.bonf.sig = p.bonf < 0.05,
         p.bonf.cat = cut(p.bonf,
                          breaks = c(1,0.05, 0.01, 0.005, 0),
                          labels = c("<= 0.005","<= 0.01", "<= 0.05", "> 0.05")
         ),
         p.fdr = p.adjust(p.value, "fdr"),
         log.p.bonf = -log10(p.bonf))

#### Cell fraction plots ####
cf_long <- pd %>% 
  select(PrimaryDx, Experiment, BrainRegion, 
         Astro, Endo, Macro, Micro, Mural, Oligo, OPC, Tcell, Excit, Inhib) %>%
  pivot_longer(!c(PrimaryDx, Experiment, BrainRegion), names_to = "CellFraction")

cf_boxplots <- ggplot(cf_long, aes(x = BrainRegion, y = value, fill = PrimaryDx, color = Experiment)) +
  geom_boxplot() +
  facet_wrap(~CellFraction, scales = "free_y", nrow = 3) +
  scale_color_manual(values = c(psychENCODE_MDD = "black", psychENCODE_BP = "skyblue")) +
  scale_fill_manual(values = mdd_Dx_colors) +
  theme_bw() +
  theme(text = element_text(size=15), legend.position="top") 

ggsave(cf_boxplots, file = here("differential_expression","plots","model_cellFraction_boxplots.png"), width = 20, height = 15)

