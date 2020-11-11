library(SummarizedExperiment)
library(jaffelab)
library(here)
library(reshape2)
library(purrr)
library(tidyverse)
library(broom)
library(sessioninfo)

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


#### Specific cell types ####
## load and melt prop data
load(here("deconvolution","data","prop_sacc.Rdata"), verbose = TRUE)
load(here("deconvolution","data","prop_amyg.Rdata"), verbose = TRUE)

prop_sacc <- melt(est_prop_sacc$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop = value)

prop_amyg <- melt(est_prop_amyg$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop = value)

## Add top 40 results
load(here("deconvolution","data","prop_top40_sacc.Rdata"), verbose = TRUE)
load(here("deconvolution","data","prop_top40_amyg.Rdata"), verbose = TRUE)

prop_sacc <- melt(est_prop_top40_sacc$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop_top40 = value) %>%
  right_join(prop_sacc,by = c("sample", "cell_type"))

prop_amyg <- melt(est_prop_top40_amyg$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop_top40 = value)%>%
  right_join(prop_amyg,by = c("sample", "cell_type"))

prop_sacc$cell_type <- factor(prop_sacc$cell_type, 
                              levels = c("Astro","Micro","Oligo","OPC",
                                         "Inhib.1","Inhib.2",
                                         "Excit.1","Excit.2","Excit.3","Excit.4"))

prop_amyg$cell_type <- factor(prop_amyg$cell_type, 
                              levels = c("Astro","Micro","Oligo","OPC",
                                         "Inhib.1","Inhib.2","Inhib.3","Inhib.4","Inhib.5",
                                         "Excit.1","Excit.2","Excit.3"))

#### broad data ####
load(here("deconvolution","data","prop_broad_sacc.Rdata"), verbose = TRUE)
load(here("deconvolution","data","prop_broad_amyg.Rdata"), verbose = TRUE)

prop_broad_sacc <- melt(est_prop_broad_sacc$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop = value) 

prop_broad_amyg <- melt(est_prop_broad_amyg$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop = value) 

load(here("deconvolution","data","prop_top40_broad_sacc.Rdata"), verbose = TRUE)
load(here("deconvolution","data","prop_top40_broad_amyg.Rdata"), verbose = TRUE)

prop_broad_sacc <- melt(est_prop_top40_broad_sacc$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop_top40 = value) %>%
  right_join(prop_broad_sacc,by = c("sample", "cell_type"))

prop_broad_amyg <- melt(est_prop_top40_broad_amyg$Est.prop.weighted) %>%
  rename(sample = Var1, cell_type = Var2, prop_top40 = value)%>%
  right_join(prop_broad_amyg,by = c("sample", "cell_type"))

prop_broad_sacc$cell_type <- factor(prop_broad_sacc$cell_type, 
                                    levels = c("Astro","Micro","Oligo","OPC","Inhib","Excit"))

prop_broad_amyg$cell_type <- factor(prop_broad_amyg$cell_type, 
                                    levels = c("Astro","Micro","Oligo","OPC","Inhib","Excit"))

## Arrange 
prop_sacc <- prop_sacc %>% arrange(cell_type)
prop_amyg <- prop_amyg %>% arrange(cell_type)
prop_broad_sacc <- prop_broad_sacc %>% arrange(cell_type)
prop_broad_amyg <- prop_broad_amyg %>% arrange(cell_type)

#### Add pd cols ####
pd_sacc <- pd_sacc %>% 
  select(PrimaryDx, AgeDeath) %>% 
  rownames_to_column("sample")

pd_amyg <- pd_amyg %>% 
  select(PrimaryDx, AgeDeath) %>% 
  rownames_to_column("sample")

## join 
prop_sacc <- prop_sacc %>%
  left_join(pd_sacc, by = "sample") 

prop_amyg <- prop_amyg %>%
  left_join(pd_amyg, by = "sample") 

prop_broad_sacc <- prop_broad_sacc %>%
  left_join(pd_sacc, by = "sample") 

prop_broad_amyg <- prop_broad_amyg %>%
  left_join(pd_amyg, by = "sample") 

save(prop_sacc, prop_amyg, prop_broad_sacc, prop_broad_amyg, 
     file = here("Deconvolution","data","prop_long.Rdata"))

# sgejobs::job_single('deconvo_data_prep', create_shell = TRUE, queue= 'bluejay', memory = '5G', command = "Rscript deconvo_data_prep.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessionsession_info()
