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
pd = colData(rse_gene)

table(pd$BrainRegion, pd$Sex, pd$PrimaryDx)
# , ,  = MDD
# 
# 
# F   M
# Amygdala  77 154
# sACC      78 150
# 
# , ,  = Control
# 
# 
# F   M
# Amygdala  38 149
# sACC      38 162
# 
# , ,  = Bipolar
# 
# 
# F   M
# Amygdala  45  77
# sACC      51  72


## Add common_feature_id, common_gene_symbol, common_gene_ID
rowData(rse_gene)$common_feature_id <- rowData(rse_gene)$gencodeID
rowData(rse_gene)$common_gene_symbol <- rowData(rse_gene)$Symbol
rowData(rse_gene)$common_gene_id <- rowData(rse_gene)$gencodeID

#### Split By Dx, Region, and Sex ####
split_samples <- function(rse){
  ## split Region
  rse_split <- map(list(amyg = "Amygdala", sacc = "sACC"), function(region){
    rse_r <- rse[,rse$BrainRegion == region]
    
    ## split Dx
    rse_d <- map(list(MDD = "MDD", BPD = "Bipolar"), function(dx){
      rse_r <- rse_r[,rse_r$PrimaryDx == dx | rse_r$PrimaryDx == "Control"]
      ## fix levels
      rse_r$PrimaryDx <- droplevels(rse_r$PrimaryDx)
      rse_r$PrimaryDx <- relevel(rse_r$PrimaryDx, "Control")
      
      ## split Sex
      rse_r <- map(list(`F` = "F", M = "M"), ~rse_r[,rse_r$Sex == .x])
      
      return(rse_r)
    })
    
    return(rse_d)
  })
  return(rse_split)
}

rse_gene_split <- split_samples(rse_gene)
map_depth(rse_gene_split, 3, ncol)

#### Define Models ####
# qSV_mat
load(here("differential_expression","data","qSV_mat.Rdata"), verbose = TRUE)

## Same as modSep from full DE, but remove Sex

modSep <- map_depth(rse_gene_split, 3, function(rse){
  qSV_split <- qSV_mat[colnames(rse),]
  mod <- model.matrix(~PrimaryDx + AgeDeath + 
                        snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
                        mitoRate + rRNA_rate +totalAssignedGene + RIN + abs(ERCCsumLogErr), 
                      data=colData(rse))
  mod <- cbind(mod, qSV_split)
  return(mod)
})

map_depth(modSep, 3, nrow)

#### RUN DE ####
source(here("differential_expression","code","run_DE.R"))

run_DE_models <- function(rse_split, model, run_voom = TRUE, csv_prefix){
    pmap(list(rse_region = rse_split, mod_region = model, region_name = names(rse_split)),
         function(rse_region, mod_region, region_name){
           
           pmap(list(rse_dx = rse_region, mod_dx = mod_region, dx_name = names(rse_region)),
                function(rse_dx, mod_dx, dx_name){
                  message("Region: ", region_name, " Dx: ", dx_name)
                  outDE<- map2(rse_dx, mod_dx, ~run_DE(rse = .x, model = .y, coef = 2, run_voom = run_voom))

                  message("Writing csv files")
                  map2(outDE, names(outDE),
                       ~write_csv(.x, file = here("differential_expression","data","DE_csv",
                                                                   paste0(csv_prefix,"_",region_name,"_",.y,".csv"))))
                  return(outDE)
                })
         })
}


message("\n####  GENE  ####")
outGene <- run_DE_models(rse_gene_split, modSep, csv_prefix = "qSVA_MDD_gene_sex")
save(outGene, file = here("differential_expression","data","qSVA_MDD_gene_DEresults.rda"))


#sgejobs::job_single('qSV_model_DE_analysis_sex', create_shell = TRUE, memory = '80G', command = "Rscript qSV_model_DE_analysis_sex.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
