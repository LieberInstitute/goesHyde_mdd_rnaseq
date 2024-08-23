# Louise Huuki-Myers Aug 22
# compile and plot all deconvolution results for MDDseq revisions

library("tidyverse")
library("here")
library("SummarizedExperiment")
library("sessioninfo")
library("writexl")
library("DeconvoBuddies")

#### Read data ####

load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
pd <- as.data.frame(colData(rse_gene))
pd2 <- pd[,c("RNum", "BrNum", "Experiment", "BrainRegion", "PrimaryDx")]

rm("rse_gene")

est_prop <- list()

## Original results - Tran Broad, Bisque
load(here("deconvolution","data","est_prop_Bisque.Rdata"),verbose = TRUE)

est_prop$tran_broad_bisque <- est_prop_bisque$Est.prop.long |>
  separate(Sample, into = c("RNum", "Experiment"), extra = "merge", remove = FALSE) |>
  left_join(pd2) |>
  filter(PrimaryDx != "Bipolar") |>
  mutate(method = "Bisque",
         refrence_dataset = "Tran_broad") |>
  arrange(RNum, cell_type, BrainRegion)

## Revis results - Tran Broad, Bisque + hspe
load(here("deconvolution","data","revis","est_prop_Bisque_fine.Rdata"),verbose = TRUE)
# est_prop_bisque

est_prop$tran_fine_bisque <- do.call("rbind", map(est_prop_bisque, "Est.prop.long")) |>
  separate(Sample, into = c("RNum", "Experiment"), extra = "merge", remove = FALSE) |>
  left_join(pd2) |>
  mutate(method = "Bisque",
         refrence_dataset = "Tran_fine") |>
  filter(PrimaryDx != "Bipolar") |>
  arrange(RNum, cell_type, BrainRegion)

## hspe
load(here("deconvolution","data","revis","est_prop_hspe_fine.Rdata"),verbose = TRUE)
# est_prop_hspe

est_prop_hspe <- map(est_prop_hspe, ~as.data.frame(.x$estimates) |>
                       rownames_to_column("Sample")|>
                       pivot_longer(!Sample, names_to = "cell_type", values_to = "prop")) 

est_prop$tran_fine_hspe <- do.call("rbind", est_prop_hspe) |>
  separate(Sample, into = c("RNum", "Experiment"), extra = "merge", remove = FALSE) |>
  left_join(pd2) |>
  mutate(method = "hspe",
         refrence_dataset = "Tran_fine") |>
  filter(PrimaryDx != "Bipolar") |>
  arrange(RNum, cell_type, BrainRegion)

## Revis results - BICCN data, hspe
est_prop_fn <- map(c(amyg = "Amygdala", sacc = "sACC"), ~here("deconvolution","BICCN_data", paste0("est_prop_BICCN_hspe_",.x,".Rdata")))
est_prop_hspe <- map(est_prop_fn, ~get(load(.x)))

est_prop$BICCN_hspe <- do.call("rbind",
                         map2(est_prop_hspe, c(amyg = "Amygdala", sacc = "sACC"), ~as.data.frame(.x$estimates) |>
                                rownames_to_column("Sample")|>
                                pivot_longer(!Sample, names_to = "cell_type", values_to = "prop") |>
                                separate(Sample, into = c("RNum", "Experiment"), extra = "merge", remove = FALSE) |>
                                left_join(pd2) |>
                                mutate(BrainRegion = .y,
                                       method = "hspe",
                                       refrence_dataset = "BICCN"))) |>
  filter(PrimaryDx != "Bipolar") |>
  arrange(RNum, cell_type, BrainRegion)

## check data
map(est_prop, ~.x |> dplyr::count(BrainRegion, cell_type) |> group_by(BrainRegion) |> dplyr::slice(1))

map(est_prop, colnames)
map(est_prop, ~.x |> select(RNum, BrNum, cell_type, method, refrence_dataset))

## write data to xlsx
write_xlsx(est_prop, here("deconvolution", "data", "MDDseq_deconvolution_est_prop.xlsx"))

#### Plot composition bar plots ####
plot_dir = here("deconvolution", "plots", "revis_ALL")


## Tran broad Bisque
ct_colors_tran_broad <- c(Astro =	"#d63468",
                         Endo =	"#d8a06d",
                         Macro = "purple",
                         Micro =	"#db7972",
                         Mural=	"#ad7e2e",
                         OPC=	"#da6f33",
                         Oligo=	"#cc3b33",
                         Tcell=	"brown4",
                         Excit =	"#abb739",
                         Inhib =	"#51aeda")

unique(est_prop$tran_broad_bisque$cell_type)
# [1] "Astro" "Endo"  "Excit" "Inhib" "Macro" "Micro" "Mural" "OPC"   "Oligo" "Tcell"
setequal(names(ct_colors_tran_broad), unique(est_prop$tran_broad_bisque$cell_type))

est_prop$tran_broad_bisque$cell_type <- factor(est_prop$tran_broad_bisque$cell_type, 
                                               levels = names(ct_colors_tran_fine))

comp_bar_tran_broad <- plot_composition_bar(
  prop_long = est_prop$tran_broad_bisque |> mutate(region_dx = paste(BrainRegion, PrimaryDx)),
  sample_col = "RNum",
  x_col = "region_dx",
  prop_col = "prop",
  ct_col = "cell_type",
  add_text = TRUE,
  min_prop_text = 0.02) +
  theme_bw() +
  labs(title = "Tran et al. - Broad Cell Types (Bisque)", x = "Brain Region + Primary Dx")

ggsave(comp_bar_tran_broad, filename = here(plot_dir, "comp_bar_tran_broad.png"))
  

## Tran fine Bisque
ct_colors_tran_fine <- c(Astro_A=	"#d63468",
                      Astro_B=	"#973a4d",
                      Micro=	"#db7972",
                      Oligo_A=	"#cc3b33",
                      Oligo_B=	"#914526",
                      OPC=	"#da6f33",
                      Endo=	"#d8a06d",
                      Mural=	"#ad7e2e",
                      Oligo=	"#dca539",
                      Tcell=	"brown4",
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

est_prop$tran_fine_bisque$cell_type <- gsub("amy_|sacc_", "", est_prop$tran_fine_bisque$cell_type)
unique(est_prop$tran_fine_bisque$cell_type)

## Tran fine hspe
est_prop$tran_fine_hspe$cell_type <- gsub("amy_|sacc_", "", est_prop$tran_fine_hspe$cell_type)
unique(est_prop$tran_fine_bisque$cell_type)

setequal(names(ct_colors_tran_fine), unique(est_prop$tran_fine_hspe$cell_type))
setequal(names(ct_colors_tran_fine), unique(est_prop$tran_fine_bisque$cell_type))

est_prop$tran_fine_bisque$cell_type <- factor(est_prop$tran_fine_bisque$cell_type, 
                                               levels = names(ct_colors_tran_fine))
est_prop$tran_fine_hspe$cell_type <- factor(est_prop$tran_fine_hspe$cell_type, 
                                               levels = names(ct_colors_tran_fine))
## BICCN 
unique(est_prop$BICCN_hspe$cell_type)

est_prop$BICCN_hspe$cell_type <- factor(est_prop$BICCN_hspe$cell_type, 
                                        levels = c("Astrocyte", "Microglia","Oligodendrocyte","Oligodendrocyte_precursor",
                                                   "Amygdala_excitatory", "CGE_interneuron", "Deep_layer_corticothalamic_and_6b",
                                                   "Deep_layer_intratelencephalic", "Eccentric_medium_spiny_neuron",
                                                   "LAMP5_LHX6_and_Chandelier","Medium_spiny_neuron",
                                                   "MGE_interneuron", "Upper_layer_intratelencephalic","Deep_layer_near_projecting"))


comp_bar_plots <- map(est_prop, ~plot_composition_bar(
  prop_long = .x |> mutate(region_dx = paste(BrainRegion,"\n", PrimaryDx)),
  sample_col = "RNum",
  x_col = "region_dx",
  prop_col = "prop",
  ct_col = "cell_type",
  add_text = TRUE,
  min_prop_text = 0.02) +
  theme_bw() +
  labs(x = "Brain Region + Primary Dx"))

walk2(comp_bar_plots, names(comp_bar_plots), ~ggsave(.x, filename = here(plot_dir, paste0("comp_bar_", .y,".png"))))

