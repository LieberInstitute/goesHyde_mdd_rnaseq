#####################################
# from /dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/R_scripts/DEG_scripts/qSV_model_analysis_DLPFC_PTSD.R
#####################################

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)
library(edgeR)
library(sessioninfo)
library(here)
library(tidyverse)

##### Load rse data, examine ####
#to see more column
options("width"=200)

#load objects
load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)
load(here('exprs_cutoff','rse_exon.Rdata'), verbose=TRUE)
load(here('exprs_cutoff','rse_jxn.Rdata'), verbose=TRUE)
load(here('exprs_cutoff','rse_tx.Rdata'), verbose=TRUE)

pd = colData(rse_gene)

table(pd$Experiment)
# psychENCODE_MDD  psychENCODE_BP 
# 588             503 

table(pd$BrainRegion)
# Amygdala     sACC 
# 540      551

table(pd$BrainRegion, pd$Sex)
# F   M
# Amygdala 160 380
# sACC     167 384


table(pd$Experiment, pd$PrimaryDx)
#                 MDD Control Bipolar
# psychENCODE_MDD 459     129       0
# psychENCODE_BP    0     258     245

table(pd$BrainRegion,pd$PrimaryDx)
#           MDD Control Bipolar
# Amygdala 231     187     122
# sACC     228     200     123


summary(pd$AgeDeath)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 17.37   34.62   47.21   46.58   55.86   95.27 

#### Split By Dx  & Region ####

split_samples <- function(rse){
  rse_split <- map(list(amyg = "Amygdala", sacc = "sACC"), function(region){
    rse_r <- rse[,rse$BrainRegion == region]
    rse_d <- map(list(MDD = "MDD", BPD = "Bipolar"), function(dx){
      rse_r <- rse_r[,rse_r$PrimaryDx == dx | rse_r$PrimaryDx == "Control"]
      ## fix levels
      rse_r$PrimaryDx <- droplevels(rse_r$PrimaryDx)
      rse_r$PrimaryDx <- relevel(rse_r$PrimaryDx, "Control")
      return(rse_r)
    })
    return(rse_d)
  })
  return(rse_split)
}

rse_gene_split <- split_samples(rse_gene)

map(rse_gene_split, ~map_int(.x, ncol))
# $amyg
# MDD BPD 
# 418 309 
# 
# $sacc
# MDD BPD 
# 428 323 

rse_exon_split <- split_samples(rse_exon)
rse_jxn_split <- split_samples(rse_jxn)
rse_tx_split <- split_samples(rse_tx)

#### Define Models ####
# qSV_mat
load(here("differential_expression","data","qSV_mat.Rdata"), verbose = TRUE)

modSep <- map(rse_gene_split, function(rse_region){
  map(rse_region, function(rse_dx){
    qSV_split <- qSV_mat[colnames(rse_dx),]
    mod <- model.matrix(~PrimaryDx + AgeDeath + Sex  + 
                    snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
                    mitoRate + rRNA_rate +totalAssignedGene + RIN + ERCCsumLogErr, 
                  data=colData(rse_dx))
    mod <- cbind(mod, qSV_split)
    return(mod)
  })
})

map(modSep, ~map(.x, head))

## Add cell fractions to model
modSep_cf <- map2(rse_gene_split, modSep, function(rse_region, mod_region){
  map2(rse_region, mod_region, ~cbind(.y, colData(.x)[,c("Astro", "Endo", "Macro", "Micro", "Mural", "Oligo", "OPC", "Tcell", "Excit")]))
})

map(modSep_cf, ~map(.x, head))

## Save Models
save(modSep, modSep_cf, file = here("differential_expression","data","differental_models.Rdata"))

#### RUN DE ####
source(here("differential_expression","code","run_DE.R"))

## Gene ##

run_DE_models <- function(rse_split, model_list = list(sep = modSep, sep_cf = modSep_cf), run_voom = TRUE, csv_prefix){
  map2(model_list, names(model_list), function(mod, mod_name){
    pmap(list(rse = rse_split, mod = mod, region_name = names(rse_split)),
         function(rse, mod, region_name){
           message("Running ", mod_name, " ", region_name, ' DE:')
           outDE<- map2(rse, mod, ~run_DE(rse = .x, model = .y, coef = 2, run_voom = run_voom))
           
           map2(outDE, names(outDE),
                ~write_csv(.x, file = here("differential_expression","data","DE_csv",
                                                            paste0(csv_prefix,"_",region_name,"_",mod_name,"_",.y,".csv"))))
           return(outDE)
         })
  })
}

outGene <- run_DE_models(rse_gene_split, csv_prefix = "qSVA_MDD_gene")
save(outGene, file = here("differential_expression","data","qSVA_MDD_gene_DEresults.rda"))

outExon <- run_DE_models(rse_exon_split, csv_prefix = "qSVA_MDD_exon")
save(outExon, file = here("differential_expression","data","qSVA_MDD_exon_DEresults.rda"))

outJxn <- run_DE_models(rse_jxn_split, csv_prefix = "qSVA_MDD_jxn")
save(outJxn, file = here("differential_expression","data","qSVA_MDD_jxn_DEresults.rda"))

outTx <- run_DE_models(rse_tx_split, run_voom = FALSE, csv_prefix = "qSVA_MDD_tx")
save(outTx, here("differential_expression","data","qSVA_MDD_tx_DEresults.rda"))

#sgejobs::job_single('qSV_model_DE_analysis', create_shell = TRUE, memory = '80G', command = "Rscript qSV_model_DE_analysis.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
