
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
                                                paste0(csv_prefix,"_",region_name,"_",dx_name,"_",.y,".csv"))))
                return(outDE)
              })
       })
}


# message("\n####  GENE  ####")
# outGene <- run_DE_models(rse_gene_split, modSep, csv_prefix = "qSVA_MDD_gene_sex")
# save(outGene, file = here("differential_expression","data","qSVA_MDD_gene_sex_DEresults.rda"))

#### Down Sample Males, rep x 1000 ####
table(rse_gene$BrainRegion, rse_gene$PrimaryDx, rse_gene$Sex)
#           F   M
# MDD     155 304
# Control  76 311
# Bipolar  96 149

down_sample <- table(rse_gene$PrimaryDx[rse_gene$Sex == "F"])
set.seed(12822)

my_downsample <- function(rse, control_n, case_n){
  
  control_i <- which(rse$PrimaryDx == "Control")
  case_i <- which(rse$PrimaryDx != "Control")
  
  control_i2 <- sample(x = control_i, size = control_n)
  case_i2 <- sample(x = case_i, size = case_n)
  
  rse <- rse[,c(control_i2, case_i2)]
  return(rse)
} 

run_DE_downsample <- function(rse_split, model, run_voom = TRUE, csv_prefix, reps = 100){
  pmap(list(rse_region = rse_split, mod_region = model, region_name = names(rse_split)),
       function(rse_region, mod_region, region_name){
         
         pmap(list(rse_dx = rse_region, mod_dx = mod_region, dx_name = names(rse_region)),
              function(rse_dx, mod_dx, dx_name){
                message("\n#########################################")
                message("Region: ", region_name, " Dx: ", dx_name)
                
                f_dx <- table(rse_dx$F$PrimaryDx == "Control")
                message("Down Sample down to:", appendLF = FALSE)
                print(table(rse_dx$F$PrimaryDx))
                
                outDE <- map(1:reps, function(r){
                  message("**** Rep: ", r, " ****")
                  rse_temp <- my_downsample(rse_dx$M, control_n = f_dx['TRUE'], case_n = f_dx['FALSE'])
                  mod_temp <- mod_dx$M[colnames(rse_temp),]
                  outDE <- run_DE(rse = rse_temp, model = mod_temp, coef = 2, run_voom = run_voom)
                  return(outDE)
                })

                message("Writing csv files")
                map2(outDE, 1:reps,
                     ~write_csv(.x, file = here("differential_expression","data","DE_csv_downsample",
                                                paste0(csv_prefix,"_",region_name,"_",dx_name,"_",.y,".csv"))))
                return(outDE)
                return(table(rse_temp$PrimaryDx))
              })
       })
}

message("Preforming M only Down sample: ")
outGene_downsample <- run_DE_downsample(rse_gene_split, modSep, csv_prefix = "qSVA_MDD_gene_sex", reps = 100)
save(outGene_downsample, file = here("differential_expression","data","qSVA_MDD_gene_M_downsample_DEresults.rda"))

#sgejobs::job_single('qSV_model_DE_analysis_sex', create_shell = TRUE, memory = '80G', command = "Rscript qSV_model_DE_analysis_sex.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
