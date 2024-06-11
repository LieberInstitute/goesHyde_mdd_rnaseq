
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
library(ggpubr)
library(rstatix)

plot_dir <- here("deconvolution", "plots", "revis")
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

load(here("deconvolution","data","revis","Tran_cell_type_prop.Rdata"),verbose = TRUE)
# cell_type_prop

## Load Bisque & hspe results
load(here("deconvolution","data","revis","est_prop_Bisque_fine.Rdata"),verbose = TRUE)
# est_prop_bisque
load(here("deconvolution","data","revis","est_prop_hspe_fine.Rdata"),verbose = TRUE)
# est_prop_hspe

prop_long_bisque <- do.call("rbind", map(est_prop_bisque, "Est.prop.long")) |>
  separate(Sample, into = c("RNum", "Experiment"), extra = "merge", remove = FALSE) |>
  left_join(pd2) |>
  mutate(region_dx = paste(BrainRegion, PrimaryDx),
         method = "Bisque")

est_prop_hspe <- map(est_prop_hspe, ~as.data.frame(.x$estimates) |>
                       rownames_to_column("Sample")|>
                       pivot_longer(!Sample, names_to = "cell_type", values_to = "prop"))

prop_long_hspe <- do.call("rbind", est_prop_hspe) |>
  separate(Sample, into = c("RNum", "Experiment"), extra = "merge", remove = FALSE) |>
  left_join(pd2) |>
  mutate(region_dx = paste(BrainRegion, PrimaryDx),
         method = "hspe")

# combined methods
sn_prop_long <- bind_rows(cell_type_prop$sacc$prop_global, cell_type_prop$amy$prop_global) |>
  rename(cell_type = cellType,
         sn_global_prop = prop)

sn_prop_long_sample <- bind_rows(cell_type_prop$sacc$sample, 
                                 cell_type_prop$amy$sample) |>
  rename(cell_type = cellType,
         sn_sample_prop = prop,
         BrNum = donor) |>
  mutate(BrNum = gsub("br", "Br", BrNum))

## all BrNum in sn data
sn_prop_long_sample |> 
  group_by(BrNum) |> 
  summarize(n = n()) |>
  mutate(in_mdd = BrNum %in% prop_long$BrNum)


# cell_types <- unique(gsub("sacc_|amy_", "", sn_prop_long$cell_type))
cell_type_levels <- c("Astro_A","Astro_B","Micro","Oligo_A","Oligo_B","OPC",
                      "Endo","Mural","Oligo","Tcell",
                      "Excit_A","Excit_B","Excit_C","Excit_D","Excit_E","Excit_F","Excit_G",
                      "Inhib_A","Inhib_B","Inhib_C","Inhib_D","Inhib_E","Inhib_F","Inhib_G","Inhib_H","Inhib_I","Inhib_J","Inhib_K")

prop_long <- bind_rows(prop_long_bisque, prop_long_hspe) |>
  mutate(region = gsub("amygdala", "amy",tolower(BrainRegion)))  |>
  filter(PrimaryDx != "Bipolar") |>
  left_join(sn_prop_long |> select(-n)) |>
  left_join(sn_prop_long_sample |> select(-n))  |>
  mutate(region_cell_type = cell_type,
         cell_type = factor(gsub("sacc_|amy_", "", cell_type),
                            levels = cell_type_levels))

prop_long |> count(cell_type)

prop_long  |>
  group_by(method, region) |>
  summarize(cor = cor(sn_global_prop, prop))
# method region   cor
# <chr>  <chr>  <dbl>
# 1 Bisque amy    0.883
# 2 Bisque sacc   0.701
# 3 hspe   amy    0.162
# 4 hspe   sacc   0.340

prop_long  |>
  filter(!is.na(sn_sample_prop)) |>
  group_by(method, region) |>
  summarize(cor = cor(sn_sample_prop, prop))
# method region   cor
# <chr>  <chr>  <dbl>
# 1 Bisque amy    0.879
# 2 Bisque sacc   0.680
# 3 hspe   amy    0.368
# 4 hspe   sacc   0.517

#### colors ####

cell_type_colors <- c(Astro_A=	"#d63468",
                      Astro_B=	"#973a4d",
                      Micro=	"#db7972",
                      Oligo_A=	"#cc3b33",
                      Oligo_B=	"#914526",
                      OPC=	"#da6f33",
                      Endo=	"#d8a06d",
                      Mural=	"#ad7e2e",
                      Oligo=	"#dca539",
                      Tcell=	"#81642b",
                      Excit_A=	"#abb739",
                      Excit_B=	"#5b6b21",
                      Excit_C=	"#92a556",
                      Excit_D=	"#5db645",
                      Excit_E=	"#35773f",
                      Excit_F=	"#58bf7e",
                      Excit_G=	"#4db5a0",
                      Inhib_A=	"#51aeda",
                      Inhib_B=	"#4771b4",
                      Inhib_C=	"#5c6ddc",
                      Inhib_D=	"#a093dd",
                      Inhib_E=	"#694a96",
                      Inhib_F=	"#803db2",
                      Inhib_G=	"#c772da",
                      Inhib_H=	"#8d3f77",
                      Inhib_I=	"#d788bc",
                      Inhib_J=	"#c83e9b",
                      Inhib_K=	"#db658f")

#### Boxplots ####
## overall
prop_boxplots <- prop_long |>
  ggplot(aes(x = cell_type, y = prop, fill = cell_type)) +
  geom_boxplot() +
  scale_fill_manual(values = cell_type_colors) +
  facet_grid(method~BrainRegion, scales = "free_x") +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(prop_boxplots, filename= here(plot_dir, "cellType_boxplot.png"))

prop_boxplots_bisque <- prop_long |>
  filter(method == "Bisque")|>
  ggplot(aes(x = cell_type, y = prop, fill = cell_type)) +
  geom_boxplot() +
  scale_fill_manual(values = cell_type_colors) +
  facet_wrap(~BrainRegion, scales = "free_x") +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(prop_boxplots_bisque, filename= here(plot_dir, "bisque_cellType_boxplot.png"), width = 10)

boxplot_dx <- prop_long |>
  ggplot(aes(x = cell_type, y = prop, fill = PrimaryDx)) +
  geom_boxplot() +
  labs(x = "Cell Type", y = "Proportion", fill ='Primary Dx') +
  theme_bw(base_size = 15)+
  scale_fill_manual(values = mdd_Dx_colors)+
  facet_grid(method~BrainRegion, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "top")

ggsave(boxplot_dx, filename = here(plot_dir,"cellType_boxplot_Dx.png"), width = 10)
# ggsave(boxplot_dx, filename = here(plot_dir,"cellType_boxplot_Dx.pdf"), width = 10)

boxplot_dx <- prop_long |>
  filter(method == "Bisque") |>
  ggplot(aes(x = cell_type, y = prop, fill = PrimaryDx)) +
  geom_boxplot() +
  labs(x = "Cell Type", y = "Proportion", fill ='Primary Dx') +
  theme_bw(base_size = 15)+
  scale_fill_manual(values = mdd_Dx_colors)+
  facet_wrap(~BrainRegion, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "top")

ggsave(boxplot_dx, filename = here(plot_dir,"bisque_cellType_boxplot_Dx.png"), width = 10)
# ggsave(boxplot_dx, filename = here(plot_dir,"cellType_boxplot_Dx.pdf"), width = 10)

## scatter plot between methods
prop_method_scatter <- prop_long |>
  select(Sample, BrainRegion, cell_type, method, prop) |>
  pivot_wider(names_from = "method", values_from = "prop") |>
  ggplot(aes(Bisque, hspe, color = cell_type)) +
  geom_point() +
  geom_abline() +
  scale_color_manual(values = cell_type_colors) +
  facet_wrap(~BrainRegion) +
  theme_bw() 
  
ggsave(prop_method_scatter, filename = here(plot_dir, "method_scatter.png"), width = 10)

#### Composition pltos ####

composition_bar_bisque <- plot_composition_bar(
  prop_long|>
    filter(method == "Bisque"), 
  sample_col = "Sample",
  x_col = "region_dx",
  min_prop_text = 0.02
) + 
  scale_fill_manual(values = cell_type_colors) +
  theme_bw()

ggsave(composition_bar_bisque, filename = here(plot_dir, "bisque_compostion_bar.png"))

## Preform t-tests
y_position <- prop_long %>%
  group_by(region_cell_type) %>%
  summarise(y.position = max(prop) + .1*max(prop))

prop_t_test_dx <- prop_long |> 
  filter(method == "Bisque") |>
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
# 1 FALSE              32
# 2 TRUE               11

prop_t_test_dx |> filter(p.bonf < 0.05) |> arrange(region_cell_type)
# region_cell_type .y.   group1 group2         p        p.adj p.format p.signif method        FDR  p.bonf p.signif.bonf
# <chr>            <chr> <chr>  <chr>      <dbl>        <dbl> <chr>    <chr>    <chr>       <dbl>   <dbl> <chr>        
#   1 amy_Astro_B      prop  MDD    Control 1.73e- 6 0.000074     1.7e-06  ****     T-test    1.24e-5 7.43e-5 ***          
#   2 amy_Micro        prop  MDD    Control 7.79e- 7 0.000033     7.8e-07  ****     T-test    6.70e-6 3.35e-5 ***          
#   3 amy_Mural        prop  MDD    Control 1.61e- 4 0.0069       0.00016  ***      T-test    8.64e-4 6.91e-3 **           
#   4 amy_Tcell        prop  MDD    Control 2.06e- 4 0.0089       0.00021  ***      T-test    9.86e-4 8.88e-3 **           
#   5 sacc_Astro_A     prop  MDD    Control 9.08e- 5 0.0039       9.1e-05  ****     T-test    5.58e-4 3.90e-3 ***          
#   6 sacc_Excit_A     prop  MDD    Control 5.30e- 4 0.023        0.00053  ***      T-test    2.07e-3 2.28e-2 *            
#   7 sacc_Excit_F     prop  MDD    Control 2.32e-10 0.00000001   2.3e-10  ****     T-test    4.98e-9 9.97e-9 ***          
#   8 sacc_Inhib_C     prop  MDD    Control 1.75e- 9 0.000000075  1.8e-09  ****     T-test    2.51e-8 7.54e-8 ***          
#   9 sacc_Inhib_H     prop  MDD    Control 2.98e- 4 0.013        0.00030  ***      T-test    1.28e-3 1.28e-2 *            
#   10 sacc_Micro       prop  MDD    Control 4.81e-11 0.0000000021 4.8e-11  ****     T-test    2.07e-9 2.07e-9 ***          
#   11 sacc_Oligo_B     prop  MDD    Control 1.53e- 8 0.00000066   1.5e-08  ****     T-test    1.65e-7 6.60e-7 ***

write.csv(prop_t_test_dx |>
            select(-y.position),
          file = here("deconvolution","data","revis","deconvolution_t_test_fine.csv"))

## TODO mutate with region
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


## Sex
prop_t_test_sex <- prop_long |> group_by(region_cell_type) |>
  do(compare_means(prop ~ Sex, data = ., method = "t.test")) |>
  ungroup() |>
  mutate(FDR = p.adjust(p, "fdr"),
         p.bonf = p.adjust(p, "bonf"),
         p.bonf.anno = paste0(round(p.bonf,3), ifelse(p.bonf < 0.05, "*",""))) |>
  left_join(y_position)

prop_t_test_sex |> count(p.bonf < 0.05)
# `p.bonf < 0.05`     n
# <lgl>           <int>
# 1 FALSE              39
# 2 TRUE                4



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
