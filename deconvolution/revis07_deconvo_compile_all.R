# Louise Huuki-Myers Aug 22
# compile and plot all deconvolution results for MDDseq revisions

library("tidyverse")
library("here")
library("SummarizedExperiment")
library("sessioninfo")
library("writexl")

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
  arrange(RNum, cell_type, BrainRegion)

map(est_prop, ~.x |> dplyr::count(BrainRegion, cell_type))

map(est_prop, colnames)

## write data
write_xlsx(est_prop, here("deconvolution", "data", "MDDseq_deconvolution_est_prop.xlsx"))

