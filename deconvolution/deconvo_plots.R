
library(SummarizedExperiment)
library(RColorBrewer)
library(tidyverse)
library(reshape2)
library(broom)
library(viridis)
library(DeconvoBuddies)
library(here)
library(sessioninfo)
library(patchwork)
library(ggpubr)
library(rstatix)

## Load colors and plotting functions
load(here("data","MDD_colors.Rdata"), verbose = TRUE)
# mdd_Dx_colors
# mdd_BrainRegion_colors
# mdd_Sex_colors
# mdd_dataset_colors

## Load Bisque results & res data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
load(here("deconvolution","data","est_prop_Bisque.Rdata"),verbose = TRUE)
load(here("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/cell_colors.Rdata"),verbose = TRUE)
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/sce_pan.v2.Rdata", verbose = TRUE)

pd <- as.data.frame(colData(rse_gene))
pd2 <- pd[,c("RNum", "BrNum", "BrainRegion","Sex", "PrimaryDx", "Experiment")]

long_prop <- est_prop_bisque$Est.prop.long %>%
  separate(Sample, into = c("RNum", "Experiment"), extra = "merge", remove = FALSE) %>%
  left_join(pd2) %>%
  mutate(region_dx = paste(BrainRegion, PrimaryDx),
         region_cell_type = paste(BrainRegion, cell_type))

long_prop$cell_type <- factor(long_prop$cell_type, levels = levels(sce_pan$cellType.Broad))

long_prop %>% filter(prop > 0.99)
# Sample                 RNum   Experiment      cell_type  prop BrNum  BrainRegion Sex   PrimaryDx
# <chr>                  <chr>  <chr>           <fct>     <dbl> <chr>  <chr>       <chr> <fct>    
# 1 R17828_psychENCODE_MDD R17828 psychENCODE_MDD Oligo      1    Br5914 sACC        F     MDD      
# 2 R15043_psychENCODE_BP  R15043 psychENCODE_BP  Oligo      1    Br5712 Amygdala    M     Control  
# 3 R15093_psychENCODE_BP  R15093 psychENCODE_BP  Inhib      1.00 Br6016 Amygdala    F     Bipolar 

## Boxplots
boxplot_dx <- long_prop %>%
  ggplot(aes(x = cell_type, y = prop, fill = PrimaryDx)) +
  geom_boxplot() +
  labs(x = "Cell Type", y = "Proportion", fill ='Primary Dx') +
  theme_bw(base_size = 15)+
  scale_fill_manual(values = mdd_Dx_colors)+
  facet_wrap(~BrainRegion) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(boxplot_dx, filename = here("deconvolution", "plots", "bisque_cellType_boxplot_Dx.png"), width = 10)
ggsave(boxplot_dx, filename = here("deconvolution", "plots", "bisque_cellType_boxplot_Dx.pdf"), width = 10)

boxplot_region <- long_prop %>%
  ggplot(aes(x = cell_type, y = prop, fill = BrainRegion)) +
  geom_boxplot() +
  labs(x = "Cell Type", y = "Proportion", fill ='Brain Region') +
  # scale_fill_manual(values = region_colors)+
  theme_bw(base_size = 15)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(boxplot_region, filename = here("deconvolution", "plots", "bisque_cellType_boxplot_region.png"))
ggsave(boxplot_region, filename = here("deconvolution", "plots", "bisque_cellType_boxplot_region.pdf"))

## Preform t-tests
y_position <- long_prop %>%
  group_by(region_cell_type) %>%
  summarise(y.position = max(prop) + .1*max(prop))

prop_t_test_dx <- long_prop %>% 
  do(compare_means(prop ~ PrimaryDx, data = ., method = "t.test", p.adjust.method = "bonferroni", group.by = "region_cell_type")) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(p, "fdr"),
         p.bonf = p.adjust(p, "bonf"),
         p.signif.bonf = case_when(p.bonf < 0.005 ~ "***",
                   p.bonf < 0.01 ~"**",
                   p.bonf < 0.05 ~"*",
                   TRUE~""),
         p.bonf.anno = paste0(round(p.bonf, 3),p.signif.bonf),
         y.adj = rep(c(0.05, 0, 0.1),20)) %>%
  left_join(y_position) %>%
  mutate(y.position = y.position - (y.adj*y.position))

prop_t_test_dx %>% count(group1, group2)

# dx_df <- data.frame(n1 = dx[c("Control", "MDD","MDD")],
#            n1 = dx[ c("Bipolar", "Bipolar", "Control")])
# colnames(dx_df) <- c("group1","n1","group2","n2")
# dx_df$df <- (dx_df$n1 + dx_df$n2 - 2)
# 
# prop_t_test_dx<- prop_t_test_dx %>% left_join(dx_df %>% select(-n1,-n2)) %>% mutate(t = qt(p, df))
# 
# summary(prop_t_test_dx$t)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -4.9128 -2.8452 -1.7731 -1.7150 -0.3765  1.9571 
t_check <- long_prop %>% 
  group_by(region_cell_type, cell_type, BrainRegion) %>% 
  summarise(mean_prop = mean(prop)) %>%
  left_join(prop_t_test_dx %>% select(region_cell_type, group1, group2, t))

t_check_scatter <- t_check %>% ggplot(aes(mean_prop, t, color = cell_type, shape = BrainRegion)) + geom_point() + facet_grid(group1~group2)
ggsave(t_check_scatter, filename = here("deconvolution", "plots","t_check_scatter.png"))

prop_t_test_dx %>% count(p.bonf < 0.05)
# `p.bonf < 0.05`     n
# <lgl>           <int>
# 1 FALSE              46
# 2 TRUE               14

prop_t_test_dx %>% filter(p.bonf < 0.05)

dx_comparisons <- list( c("MDD", "Control"), c("Bipolar", "Control"), c("MDD", "Bipolar") )

ggbox <- ggboxplot(long_prop, x = "PrimaryDx", y = "prop", fill = "PrimaryDx", facet.by = "region_cell_type", scales = "free_y") +
  scale_fill_manual(values = mdd_Dx_colors)+
  # stat_compare_means(aes(label=..p.adj..),
  #                    comparisons = dx_comparisons, method = "t.test", 
  #                    p.adjust.method = "bonferroni", 
  #                    group.by = "region_cell_type")+
  stat_pvalue_manual(prop_t_test_dx, label = "p.bonf.anno", color = "p.signif.bonf")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(ggbox, filename = here("deconvolution", "plots", "ggbox_t-test_Dx.png"), width = 16, height = 16)

## Sex
prop_t_test_sex <- long_prop %>% group_by(region_cell_type) %>%
  do(compare_means(prop ~ Sex, data = ., method = "t.test")) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(p, "fdr"),
         p.bonf = p.adjust(p, "bonf"),
         p.bonf.anno = paste0(round(p.bonf,3), ifelse(p.bonf < 0.05, "*",""))) %>%
  left_join(y_position)

prop_t_test_sex %>% count(p.bonf < 0.05)
# `p.bonf < 0.05`     n
# <lgl>           <int>
#   1 FALSE              20

## Manual
ggbox_sex <- ggboxplot(long_prop, x = "Sex", y = "prop", fill = "Sex", facet.by = "region_cell_type", scales = "free_y") + 
  stat_pvalue_manual(prop_t_test_sex, label = "p.bonf.anno", color = "blue")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(ggbox_sex, filename = here("deconvolution", "plots", "ggbox_t-test_Sex.png"), width = 15, height = 15)


## Check out effect size and power
## Amyg
d = 0.2

prop_sd <- long_prop %>% group_by(BrainRegion, cell_type) %>% summarise(sd = sd(prop))

effect_size <- long_prop %>% group_by(BrainRegion, Sex, cell_type) %>% 
  summarise(mean = mean(prop)) %>%
  pivot_wider(names_from = Sex, values_from = mean) %>%
  left_join(prop_sd) %>%
  mutate(d = abs(`F` - M)/sd)

effect_size %>%
  group_by(BrainRegion) %>%
  summarise(min(d), max(d))

# BrainRegion `min(d)` `max(d)`
# <chr>          <dbl>    <dbl>
# 1 Amygdala      0.0313    0.254
# 2 sACC          0.0132    0.254

pwr.t2n.test(n1 = 160, n2= 380, d = 0.0313 ,sig.level = 0.05)
# power = 0.06268394
pwr.t2n.test(n1 = 160, n2= 380, d = 0.254 ,sig.level = 0.05)
# power = 0.7674
## sACC 
pwr.t2n.test(n1 = 167, n2= 384, d = 0.0132 ,sig.level = 0.05)
# power = 0.05231809
pwr.t2n.test(n1 = 167, n2= 384, d = 0.254 ,sig.level = 0.05)
# power = 0.7809604

#### composition barplot ####
comp_barplot <- plot_composition_bar(long_prop, x_col = "BrainRegion", min_prop_text = 0.01) +
  scale_fill_manual(values = cell_colors)+
  theme_bw(base_size = 15)

ggsave(comp_barplot, filename = here("deconvolution","plots","bisque_composition_barplot.png"))
ggsave(comp_barplot, filename = here("deconvolution","plots","bisque_composition_barplot.pdf"))

comp_barplot_dx <- plot_composition_bar(long_prop, x_col = "region_dx", min_prop_text = 0.01) +
  scale_fill_manual(values = cell_colors)+
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(comp_barplot_dx, filename = here("deconvolution","plots","bisque_composition_barplot_dx.png"), height = 8)
ggsave(comp_barplot, filename = here("deconvolution","plots","bisque_composition_barplot.pdf"))


comp_barplot_sample <- plot_composition_bar(long_prop, x_col = "RNum",
                                            add_text = FALSE) +
  scale_fill_manual(values = cell_colors)+
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(comp_barplot_sample, filename = here("deconvolution","plots","bisque_composition_barplot_sample.png"), height = 8, width = 15)

#### SCE composition ####
sce_pd<- as.data.frame(colData(sce_pan)) %>%
  mutate(Sample = paste(protocol, region, donor))

cell_counts <- sce_pd  %>%
  group_by(Sample, cellType.Broad) %>%
  summarise(n_cellType = n())

cell_prop <- sce_pd  %>%
  group_by(Sample,protocol, region, donor) %>%
  summarise(total_nuc = n()) %>%
  left_join(cell_counts) %>%
  mutate(prop = n_cellType/total_nuc)

levels(cell_prop$cellType.Broad)
levels(long_prop$cell_type)

walk(c("ALL","Sample","region"), function(x){
  comp_prop_plot <- plot_composition_bar(cell_prop,
                                         sample_col = "Sample",
                                         x_col = x,
                                         ct_col = "cellType.Broad") +
    scale_fill_manual(values = cell_colors)
  ggsave(comp_prop_plot, filename = here("deconvolution", "plots",paste0("sce_",x,"_prop.png")))
})


#### Cor with qSV ####
load(here("differential_expression", "data","qSV_mat.Rdata"), verbose = TRUE)
dim(qSV_mat)

qSV_long <- melt(qSV_mat) %>% rename(Sample = Var1, qSV = Var2, qSV_value = value)

## Bind with qSV table
est_prop_qsv <- left_join(long_prop, qSV_long, by = "Sample")
levels(est_prop_qsv$cell_type)

#### Calculate p-values ####
prop_qSV_fit <- est_prop_qsv %>% group_by(cell_type, qSV, BrainRegion) %>%
  do(fitQSV = tidy(lm(prop ~ qSV_value, data = .))) %>%
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

levels(prop_qSV_fit$p.bonf.cat)
prop_qSV_fit %>% count(BrainRegion, p.bonf.cat)
# BrainRegion p.bonf.cat     n
# <chr>       <fct>      <int>
# 1 Amygdala    <= 0.005      53
# 2 Amygdala    <= 0.01        3
# 3 Amygdala    <= 0.05       12
# 4 Amygdala    > 0.05       192
# 5 sACC        <= 0.005      58
# 6 sACC        <= 0.01        4
# 7 sACC        <= 0.05        9
# 8 sACC        > 0.05       189

#### Tile plots ####
my_breaks <- c(0.05, 0.01, 0.005, 0)

sig_colors <- c(rev(viridis_pal(option = "magma")(4)))
names(sig_colors) <- levels(prop_qSV_fit$p.bonf.cat)

tile_plot_val <- prop_qSV_fit %>%
  ggplot(aes(x = qSV, y = cell_type, fill = log.p.bonf)) +
  geom_tile(color = "grey") +
  geom_text(aes(label = ifelse(p.bonf.sig, format(round(-log10(p.bonf),1), nsmall = 1),""),
                color = p.bonf.cat), size = 3, fontface = "bold",
            show.legend = F)+
  scale_color_manual(values = sig_colors) +
  scale_fill_viridis(name = "-log10(p-value Bonf)", option = "magma", direction = -1) +
  facet_wrap(~BrainRegion, ncol = 1)+
  scale_y_discrete(limits = rev) +
  labs(title ="p-values cell-type prop~qSV", x = 'Cell Type', color = "p-value Bonf\nsignificance") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(plot = tile_plot_val,
       filename = here("deconvolution","plots","qSV_prop_fit_tileVal.png"),
       width = 15)

ggsave(plot = tile_plot_val,
       filename = here("deconvolution","plots","qSV_prop_fit_tileVal.pdf"),
       width = 10)

#### Create scatter plots ####
# sig_colors2 <- c(brewer.pal(3, "Set1"),"black")
sig_colors2 <- c("#440154","#31688E", "#35B779","black")
names(sig_colors2) <- levels(prop_qSV_fit$p.bonf.cat)

regions <- list(sacc = "sACC", amyg = "Amygdala")

est_prop_qsv_fit <- left_join(est_prop_qsv, prop_qSV_fit)

scatter_plot <- map(regions, ~  est_prop_qsv_fit %>%
                      filter(BrainRegion == .x) %>%
                      ggplot(aes(x = qSV_value, y = prop, color = p.bonf.cat))+
                      geom_point(size = .4, alpha = .2) +
                      facet_grid(cell_type~qSV, scales = "free")+
                      theme_bw(base_size = 10)+
                      scale_color_manual(values = sig_colors2) +
                      theme(legend.text = element_text(size = 15)) +
                      guides(color = guide_legend(override.aes = list(size=5))) +
                      labs(title = .x)
)


walk2(scatter_plot, names(scatter_plot), ~ggsave(filename = here("deconvolution","plots", paste0("qSV_cellType_scatter-",.y,".png")),
       plot = .x, width = 26, height = 10))

#### Cor w/ Covariates ####
covariate_terms <- c("rRNA_rate","RIN","overallMapRate","mitoRate","AgeDeath")
qc <- pd[,c("RNum",covariate_terms)]

qc_long <- qc %>% pivot_longer(!RNum, names_to = "Covariate")
head(qc_long)

## Bind with qc table
est_prop_qc <- left_join(long_prop, qc_long, by = "RNum")

#### Calculate p-values ####
prop_qc_fit <- est_prop_qc %>% group_by(cell_type, Covariate, BrainRegion) %>%
  do(fitQSV = tidy(lm(prop ~ value, data = .))) %>%
  unnest(fitQSV) %>%
  filter(term == "value") %>%
  mutate(p.bonf = p.adjust(p.value, "bonf"),
         p.bonf.sig = p.bonf < 0.05,
         p.bonf.cat = cut(p.bonf,
                          breaks = c(1,0.05, 0.01, 0.005, 0),
                          labels = c("<= 0.005","<= 0.01", "<= 0.05", "> 0.05")
         ),
         p.fdr = p.adjust(p.value, "fdr"),
         log.p.bonf = -log10(p.bonf))

levels(prop_qc_fit$p.bonf.cat)
prop_qc_fit %>% count(BrainRegion, p.bonf.cat)
# BrainRegion p.bonf.cat     n
# <chr>       <fct>      <int>
# 1 Amygdala    <= 0.005      10
# 2 Amygdala    <= 0.05        4
# 3 Amygdala    > 0.05        36
# 4 sACC        <= 0.005      11
# 5 sACC        <= 0.01        1
# 6 sACC        <= 0.05        2
# 7 sACC        > 0.05        36

## covariate tile plots
tile_plot_val_qc <- prop_qc_fit %>%
  ggplot(aes(x = Covariate, y = cell_type, fill = log.p.bonf)) +
  geom_tile(color = "grey") +
  geom_text(aes(label = ifelse(p.bonf.sig, format(round(-log10(p.bonf),1), nsmall = 1),""),
                color = p.bonf.cat), size = 3, fontface = "bold",
            show.legend = F)+
  scale_color_manual(values = sig_colors) +
  scale_fill_viridis(name = "-log10(p-value Bonf)", option = "magma", direction = -1) +
  facet_wrap(~BrainRegion, ncol = 1)+
  scale_y_discrete(limits = rev) +
  labs(title ="p-values cell-type prop~qSV", x = 'Cell Type', color = "p-value Bonf\nsignificance") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(plot = tile_plot_val_qc,
       filename = here("deconvolution","plots","covar_prop_fit_tileVal.png"))

## Scatter plots
regions <- list(sacc = "sACC", amyg = "Amygdala")

est_prop_qc_fit <- left_join(est_prop_qc, prop_qc_fit)

scatter_plot_qc <- map(regions, ~  est_prop_qc_fit %>%
                      filter(BrainRegion == .x) %>%
                      ggplot(aes(x = value, y = prop, color = p.bonf.cat))+
                      geom_point(size = .4, alpha = .2) +
                      facet_grid(cell_type~Covariate, scales = "free")+
                      theme_bw(base_size = 10)+
                      scale_color_manual(values = sig_colors2) +
                      theme(legend.text = element_text(size = 15)) +
                      guides(color = guide_legend(override.aes = list(size=5))) +
                      labs(title = .x)
)


walk2(scatter_plot_qc, names(scatter_plot_qc), ~ggsave(filename = here("deconvolution","plots", paste0("covar_cellType_scatter-",.y,".png")),
                                                 plot = .x, width = 10, height = 10))

# sgejobs::job_single('deconvo_plots', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript deconvo_plots.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.0 Patched (2021-05-18 r80330)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2021-06-21
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date       lib source
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# backports              1.2.1    2020-12-09 [2] CRAN (R 4.1.0)
# beachmat               2.8.0    2021-05-19 [2] Bioconductor
# Biobase              * 2.52.0   2021-05-19 [2] Bioconductor
# BiocGenerics         * 0.38.0   2021-05-19 [2] Bioconductor
# BiocNeighbors          1.10.0   2021-05-19 [1] Bioconductor
# BiocParallel           1.26.0   2021-05-19 [2] Bioconductor
# BiocSingular           1.8.1    2021-06-08 [1] Bioconductor
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# bluster                1.2.1    2021-05-27 [1] Bioconductor
# broom                * 0.7.7    2021-06-13 [2] CRAN (R 4.1.0)
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.1.0)
# cli                    2.5.0    2021-04-26 [2] CRAN (R 4.1.0)
# cluster                2.1.2    2021-04-17 [3] CRAN (R 4.1.0)
# codetools              0.2-18   2020-11-04 [2] CRAN (R 4.1.0)
# colorout             * 1.2-2    2021-05-27 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.0-1    2021-05-04 [2] CRAN (R 4.1.0)
# crayon                 1.4.1    2021-02-08 [2] CRAN (R 4.1.0)
# DBI                    1.1.1    2021-01-15 [2] CRAN (R 4.1.0)
# dbplyr                 2.1.1    2021-04-06 [2] CRAN (R 4.1.0)
# DeconvoBuddies       * 0.99.0   2021-06-14 [1] Github (lahuuki/DeconvoBuddies@d889340)
# DelayedArray           0.18.0   2021-05-19 [2] Bioconductor
# DelayedMatrixStats     1.14.0   2021-05-19 [2] Bioconductor
# digest                 0.6.27   2020-10-24 [2] CRAN (R 4.1.0)
# dplyr                * 1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
# dqrng                  0.3.0    2021-05-01 [1] CRAN (R 4.1.0)
# edgeR                  3.34.0   2021-05-19 [2] Bioconductor
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  0.5.0    2021-05-25 [2] CRAN (R 4.1.0)
# farver                 2.1.0    2021-02-28 [2] CRAN (R 4.1.0)
# forcats              * 0.5.1    2021-01-27 [2] CRAN (R 4.1.0)
# fs                     1.5.0    2020-07-31 [2] CRAN (R 4.1.0)
# generics               0.1.0    2020-10-31 [2] CRAN (R 4.1.0)
# GenomeInfoDb         * 1.28.0   2021-05-19 [2] Bioconductor
# GenomeInfoDbData       1.2.6    2021-05-11 [2] Bioconductor
# GenomicRanges        * 1.44.0   2021-05-19 [2] Bioconductor
# ggplot2              * 3.3.4    2021-06-16 [2] CRAN (R 4.1.0)
# glue                   1.4.2    2020-08-27 [2] CRAN (R 4.1.0)
# gridExtra              2.3      2017-09-09 [2] CRAN (R 4.1.0)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# haven                  2.4.1    2021-04-23 [2] CRAN (R 4.1.0)
# here                 * 1.0.1    2020-12-13 [1] CRAN (R 4.1.0)
# hms                    1.1.0    2021-05-17 [2] CRAN (R 4.1.0)
# httr                   1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
# igraph                 1.2.6    2020-10-06 [2] CRAN (R 4.1.0)
# IRanges              * 2.26.0   2021-05-19 [2] Bioconductor
# irlba                  2.3.3    2019-02-05 [2] CRAN (R 4.1.0)
# jsonlite               1.7.2    2020-12-09 [2] CRAN (R 4.1.0)
# labeling               0.4.2    2020-10-20 [2] CRAN (R 4.1.0)
# lattice                0.20-44  2021-05-02 [3] CRAN (R 4.1.0)
# lifecycle              1.0.0    2021-02-15 [2] CRAN (R 4.1.0)
# limma                  3.48.0   2021-05-19 [2] Bioconductor
# locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
# lubridate              1.7.10   2021-02-26 [2] CRAN (R 4.1.0)
# magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
# Matrix                 1.3-4    2021-06-01 [3] CRAN (R 4.1.0)
# MatrixGenerics       * 1.4.0    2021-05-19 [2] Bioconductor
# matrixStats          * 0.59.0   2021-06-01 [2] CRAN (R 4.1.0)
# metapod                1.0.0    2021-05-19 [1] Bioconductor
# modelr                 0.1.8    2020-05-19 [2] CRAN (R 4.1.0)
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# patchwork            * 1.1.1    2020-12-17 [1] CRAN (R 4.1.0)
# pillar                 1.6.1    2021-05-16 [2] CRAN (R 4.1.0)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# pryr                   0.1.4    2018-02-18 [2] CRAN (R 4.1.0)
# ps                     1.6.0    2021-02-28 [2] CRAN (R 4.1.0)
# purrr                * 0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R6                     2.5.0    2020-10-28 [2] CRAN (R 4.1.0)
# rafalib                1.0.0    2015-08-09 [1] CRAN (R 4.1.0)
# ragg                   1.1.3    2021-06-09 [1] CRAN (R 4.1.0)
# RColorBrewer         * 1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
# Rcpp                   1.0.6    2021-01-15 [2] CRAN (R 4.1.0)
# RCurl                  1.98-1.3 2021-03-16 [2] CRAN (R 4.1.0)
# readr                * 1.4.0    2020-10-05 [2] CRAN (R 4.1.0)
# readxl                 1.3.1    2019-03-13 [2] CRAN (R 4.1.0)
# reprex                 2.0.0    2021-04-02 [2] CRAN (R 4.1.0)
# rlang                * 0.4.11   2021-04-30 [2] CRAN (R 4.1.0)
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
# rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.1.0)
# rsvd                   1.0.5    2021-04-16 [1] CRAN (R 4.1.0)
# rvest                  1.0.0    2021-03-09 [2] CRAN (R 4.1.0)
# S4Vectors            * 0.30.0   2021-05-19 [2] Bioconductor
# ScaledMatrix           1.0.0    2021-05-19 [1] Bioconductor
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
# scran                  1.20.1   2021-05-24 [1] Bioconductor
# scuttle                1.2.0    2021-05-19 [1] Bioconductor
# sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.1.0)
# sgejobs                0.99.1   2021-05-27 [1] Github (LieberInstitute/sgejobs@f5ab0ca)
# SingleCellExperiment   1.14.1   2021-05-21 [2] Bioconductor
# sparseMatrixStats      1.4.0    2021-05-19 [2] Bioconductor
# statmod                1.4.36   2021-05-10 [2] CRAN (R 4.1.0)
# stringi                1.6.2    2021-05-17 [2] CRAN (R 4.1.0)
# stringr              * 1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment * 1.22.0   2021-05-19 [2] Bioconductor
# systemfonts            1.0.2    2021-05-11 [2] CRAN (R 4.1.0)
# textshaping            0.3.5    2021-06-09 [1] CRAN (R 4.1.0)
# tibble               * 3.1.2    2021-05-16 [2] CRAN (R 4.1.0)
# tidyr                * 1.1.3    2021-03-03 [2] CRAN (R 4.1.0)
# tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
# tidyverse            * 1.3.1    2021-04-15 [2] CRAN (R 4.1.0)
# utf8                   1.2.1    2021-03-12 [2] CRAN (R 4.1.0)
# vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
# viridis              * 0.6.1    2021-05-11 [2] CRAN (R 4.1.0)
# viridisLite          * 0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
# withr                  2.4.2    2021-04-18 [2] CRAN (R 4.1.0)
# xml2                   1.3.2    2020-04-23 [2] CRAN (R 4.1.0)
# XVector                0.32.0   2021-05-19 [2] Bioconductor
# zlibbioc               1.38.0   2021-05-19 [2] Bioconductor
#
# [1] /users/lhuuki/R/4.1
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/library
