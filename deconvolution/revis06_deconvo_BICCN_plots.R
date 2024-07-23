
library("SummarizedExperiment")
library("RColorBrewer")
library("tidyverse")
library("reshape2")
library("broom")
library("viridis")
library("DeconvoBuddies")
library("here")
library("sessioninfo")
library("patchwork")
library("ggpubr")
library("rstatix")

plot_dir <- here("deconvolution", "plots", "revis_BICCN")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

## Load colors and plotting functions
load(here("data","MDD_colors.Rdata"), verbose = TRUE)
# mdd_Dx_colors
# mdd_BrainRegion_colors
# mdd_Sex_colors
# mdd_dataset_colors

## Load res and cell prop data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
pd <- as.data.frame(colData(rse_gene))
pd2 <- pd[,c("RNum", "BrNum", "BrainRegion","Sex", "PrimaryDx", "Experiment")]

## Load Bisque & hspe results
load(here("deconvolution","BICCN_data", "est_prop_BICCN_hspe_Amygdala.Rdata"), verbose = TRUE)
head(est_prop_hspe$estimates)

est_prop_fn <- map(c(amyg = "Amygdala", sacc = "sACC"), ~here("deconvolution","BICCN_data", paste0("est_prop_BICCN_hspe_",.x,".Rdata")))
map(est_prop_fn, file.exists)
est_prop_hspe <- map(est_prop_fn, ~get(load(.x)))

est_prop_hspe <- map2(est_prop_hspe, c(amyg = "Amygdala", sacc = "sACC"), ~as.data.frame(.x$estimates) |>
                       rownames_to_column("Sample")|>
                       pivot_longer(!Sample, names_to = "cell_type", values_to = "prop") |>
                        mutate(BrainRegion = .y))

cell_type_levels <- c("Astrocyte", "Microglia","Oligodendrocyte","Oligodendrocyte_precursor",
                      "Amygdala_excitatory", "CGE_interneuron", "Deep_layer_corticothalamic_and_6b",
                      "Deep_layer_intratelencephalic", "Eccentric_medium_spiny_neuron",
                      "LAMP5_LHX6_and_Chandelier","Medium_spiny_neuron",
                      "MGE_interneuron", "Upper_layer_intratelencephalic","Deep_layer_near_projecting")


prop_long <- do.call("rbind", est_prop_hspe) |>
  separate(Sample, into = c("RNum", "Experiment"), extra = "merge", remove = FALSE) |>
  left_join(pd2) |>
  mutate(region_dx = paste(BrainRegion, PrimaryDx),
         region_cell_type = paste(BrainRegion, cell_type),
         method = "hspe",
         cell_type = factor(cell_type, levels = cell_type_levels)) |>
  filter(PrimaryDx != "Bipolar")

prop_long |> count(BrainRegion, cell_type) |> print(n=25)

#### colors ####

# cell_type_colors <- c()

#### Boxplots ####
## overall
prop_boxplots <- prop_long |>
  ggplot(aes(x = cell_type, y = prop, fill = cell_type)) +
  geom_boxplot() +
  # scale_fill_manual(values = cell_type_colors) +
  facet_grid(~BrainRegion, scales = "free_x") +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(prop_boxplots, filename= here(plot_dir, "cellType_boxplot.png"))

boxplot_dx <- prop_long |>
  ggplot(aes(x = cell_type, y = prop, fill = PrimaryDx)) +
  geom_boxplot() +
  labs(x = "Cell Type", y = "Proportion", fill ='Primary Dx') +
  theme_bw(base_size = 15)+
  scale_fill_manual(values = mdd_Dx_colors)+
  facet_grid(~BrainRegion, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "top")

ggsave(boxplot_dx, filename = here(plot_dir,"cellType_boxplot_Dx.png"), width = 10)
# ggsave(boxplot_dx, filename = here(plot_dir,"cellType_boxplot_Dx.pdf"), width = 10)


#### Composition plots ####

composition_bar_Dx <- plot_composition_bar(
  prop_long, 
  sample_col = "Sample",
  x_col = "region_dx",
  min_prop_text = 0.02
) + 
  # scale_fill_manual(values = cell_type_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(composition_bar_Dx, filename = here(plot_dir, "BICCN_est_prop_compostion_bar.png"))

## Preform t-tests
y_position <- prop_long %>%
  group_by(region_cell_type) %>%
  summarise(y.position = max(prop) + .1*max(prop))

prop_t_test_dx <- prop_long |> 
  do(compare_means(prop ~ PrimaryDx, data = ., method = "t.test", p.adjust.method = "bonferroni", group.by = "region_cell_type")) |>
  ungroup() |>
  mutate(FDR = p.adjust(p, "fdr"),
         p.bonf = p.adjust(p, "bonf"),
         p.signif.bonf = case_when(p.bonf < 0.005 ~ "***",
                   p.bonf < 0.01 ~"**",
                   p.bonf < 0.05 ~"*",
                   TRUE~""),
         p.bonf.anno = paste0(round(p.bonf, 3),p.signif.bonf),
         cell_type = factor(gsub("sacc_|amy_", "", region_cell_type),
                            levels = cell_type_levels)
         ) |>
  left_join(y_position)

prop_t_test_dx |> count(p.bonf < 0.05)
# `p.bonf < 0.05`     n
# <lgl>           <int>
# 1 FALSE              21
# 2 TRUE                4

prop_t_test_dx |> filter(p.bonf < 0.05) |> arrange(region_cell_type)
#   region_cell_type            .y.   group1 group2       p   p.adj p.format p.signif method     FDR  p.bonf p.signif.bonf
#   <chr>                       <chr> <chr>  <chr>    <dbl>   <dbl> <chr>    <chr>    <chr>    <dbl>   <dbl> <chr>        
#   1 Amygdala Microglia          prop  MDD    Contr… 7.15e-7 1.80e-5 7.1e-07  ****     T-test 8.93e-6 1.79e-5 ***          
#   2 Amygdala Oligodendrocyte_p… prop  MDD    Contr… 1.93e-4 4.8 e-3 0.00019  ***      T-test 1.21e-3 4.83e-3 ***          
#   3 sACC LAMP5_LHX6_and_Chande… prop  MDD    Contr… 1.43e-6 3.60e-5 1.4e-06  ****     T-test 1.19e-5 3.58e-5 ***          
#   4 sACC Microglia              prop  MDD    Contr… 1.01e-7 2.5 e-6 1.0e-07  ****     T-test 2.52e-6 2.52e-6 *** 

write.csv(prop_t_test_dx |>
            select(-y.position),
          file = here("deconvolution","BICCN_data","deconvolution_t_test_BICCN.csv"))

prop_t_test_signif <- prop_t_test_dx |> filter(p.bonf < 0.05) |> mutate()
prop_t_test_signif |> count(region_cell_type)

ggbox <- ggpubr::ggboxplot(prop_long,
                           x = "PrimaryDx", 
                           y = "prop", 
                           fill = "PrimaryDx", 
                           facet.by = "region_cell_type", 
                           scales = "free_y") +
  scale_fill_manual(values = mdd_Dx_colors)+
  # stat_compare_means(aes(label=..p.adj..),
  #                    comparisons = dx_comparisons, method = "t.test", 
  #                    p.adjust.method = "bonferroni", 
  #                    group.by = "region_cell_type")+
  stat_pvalue_manual(prop_t_test_dx, label = "p.bonf.anno", color = "p.signif.bonf")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(ggbox, filename = here(plot_dir,  "ggbox_t-test_Dx.png"), width = 16, height = 16)


#### Region specific significant plots ####
## amy
prop_long_s_amyg <- prop_long |> 
  filter(region_cell_type %in% (prop_t_test_signif |> pull(region_cell_type)),
         BrainRegion == "Amygdala")

ggbox_s_amyg <- ggboxplot(prop_long_s_amyg,
                          x = "PrimaryDx",
                          y = "prop",
                          fill = "PrimaryDx", 
                          facet.by = "cell_type", 
                          scales = "free_y", nrow = 1) +
  scale_fill_manual(values = mdd_Dx_colors)+
  stat_pvalue_manual(prop_t_test_signif |> 
                       filter(grepl("amy_", region_cell_type)),
                     label = "p.bonf.anno", 
                     color = "p.signif.bonf")+
  theme(legend.position = "none") +
  labs(title = "Amygdala", y = "Cell Type Proportion")

ggsave(ggbox_s_amyg, filename = here(plot_dir,  "ggbox_t-test_Dx_s_Amyg.png"), width = 9, height = 5)

## sacc
prop_long_s_sacc <- prop_long |> 
  filter(region_cell_type %in% (prop_t_test_signif |> pull(region_cell_type)),
         BrainRegion == "sACC")

ggbox_s_sacc <- ggboxplot(prop_long_s_sacc,
                          x = "PrimaryDx",
                          y = "prop",
                          fill = "PrimaryDx", 
                          facet.by = "cell_type", 
                          scales = "free_y", nrow = 1) +
  scale_fill_manual(values = mdd_Dx_colors)+
  stat_pvalue_manual(prop_t_test_signif |> 
                       filter(grepl("sacc_", region_cell_type)),
                     label = "p.bonf.anno", 
                     color = "p.signif.bonf")+
  theme(legend.position = "none") +
  labs(title = "sACC", y = "Cell Type Proportion")


ggsave(ggbox_s_sacc, filename = here(plot_dir,  "ggbox_t-test_Dx_s_sacc.png"), width = 12, height = 6)

#### Input data proportion plots ####

biccn_ct_prop <- read.csv(here("deconvolution", "BICCN_data", "cell_type_prop_Amygdala.csv")) |>
  mutate(BrainRegion = "Amygdala") |>
  bind_rows(read.csv(here("deconvolution", "BICCN_data", "cell_type_prop_sACC.csv")) |>
              mutate(BrainRegion = "sACC") ) |>
  mutate(cell_type = factor(supercluster, levels = cell_type_levels))


# composition_bar_input <- plot_composition_bar(
#   biccn_ct_prop, 
#   sample_col = "donor_id",
#   x_col = "BrainRegion",
#   ct_col = "supercluster",
#   min_prop_text = 0.02
# ) + 
#   # scale_fill_manual(values = cell_type_colors) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust=1))

composition_bar_input_sample <-  biccn_ct_prop |>
  ggplot(aes(x = donor_id, y = prop, fill = cell_type)) +
  geom_col() +
  geom_text(aes(label = ifelse(prop > 0.02,format(round(prop, 3), 3),"")),
                position = ggplot2::position_stack(vjust = 0.5)) +
  facet_wrap(~BrainRegion) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(composition_bar_input_sample, filename = here(plot_dir, "BICCN_input_sample_compostion_bar.png"))


biccn_ct_prop_region <- biccn_ct_prop |>
  group_by(BrainRegion, cell_type) |>
  summarise(mean_prop = mean(prop),
            n= sum(n)) |>
  group_by(BrainRegion) |>
  mutate(prop = n/sum(n))

composition_bar_input<-  biccn_ct_prop_region |>
  ggplot(aes(x = BrainRegion, y = prop, fill = cell_type)) +
  geom_col() +
  geom_text(aes(label = ifelse(prop > 0.02,format(round(prop, 3), 3),"")),
            position = ggplot2::position_stack(vjust = 0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(composition_bar_input, filename = here(plot_dir, "BICCN_input_compostion_bar.png"))



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
