#####################################
# from /dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/R_scripts/DEG_scripts/qSV_model_analysis_DLPFC_PTSD.R
#####################################

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)
library(edgeR)
library(purrr)
library(sessioninfo)
library(here)

source("run_DE.R")

##### Load Data ####

## load rse
load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)
pd = colData(rse_gene)

## qSV data
load(here("differential_expression","data","qSV_mat.Rdata"), verbose = TRUE)

## load ilr data
load(here("deconvolution","data","est_prop_ilr.Rdata"), verbose = TRUE)
load(here("deconvolution","data","est_prop_top5.Rdata"), verbose = TRUE)
est_prop <- map(est_prop, "Est.prop.weighted")

#### Define models ####
regions <- list(sacc = "sACC", amyg = "Amygdala")

modSep <- map(regions, ~cbind(model.matrix(~PrimaryDx + AgeDeath + Sex  + mitoRate + 
                                           rRNA_rate +totalAssignedGene + RIN + ERCCsumLogErr +
                                           snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5,
                                           data= pd[pd$BrainRegion == .x,]),
                              qSV_mat[pd$BrainRegion == .x,]))

map(modSep, colnames)
# $sacc
# [1] "(Intercept)"       "PrimaryDxControl"  "PrimaryDxBipolar"  "AgeDeath"          "SexM"             
# [6] "mitoRate"          "rRNA_rate"         "totalAssignedGene" "RIN"               "ERCCsumLogErr"    
# [11] "snpPC1"            "snpPC2"            "snpPC3"            "snpPC4"            "snpPC5"           
# [16] "PC1"               "PC2"               "PC3"               "PC4"               "PC5"              
# [21] "PC6"               "PC7"               "PC8"               "PC9"               "PC10"             
# [26] "PC11"              "PC12"              "PC13"              "PC14"              "PC15"             
# [31] "PC16"              "PC17"              "PC18"              "PC19"              "PC20"             
# [36] "PC21"              "PC22"              "PC23"              "PC24"              "PC25"             
# [41] "PC26"             
# 
# $amyg
# [1] "(Intercept)"       "PrimaryDxControl"  "PrimaryDxBipolar"  "AgeDeath"          "SexM"             
# [6] "mitoRate"          "rRNA_rate"         "totalAssignedGene" "RIN"               "ERCCsumLogErr"    
# [11] "snpPC1"            "snpPC2"            "snpPC3"            "snpPC4"            "snpPC5"           
# [16] "PC1"               "PC2"               "PC3"               "PC4"               "PC5"              
# [21] "PC6"               "PC7"               "PC8"               "PC9"               "PC10"             
# [26] "PC11"              "PC12"              "PC13"              "PC14"              "PC15"             
# [31] "PC16"              "PC17"              "PC18"              "PC19"              "PC20"             
# [36] "PC21"              "PC22"              "PC23"              "PC24"              "PC25"             
# [41] "PC26"  

## Add ilr terms
modSep_ilr <- map2(modSep, est_prop_ilr[c("sacc_specific","amyg_specific")], ~cbind(.x, .y))
map(modSep_ilr, colnames)

## Add prop terms
map(est_prop, colnames)

modSep_prop <- map2(modSep, est_prop[c("sacc_specific","amyg_specific")], ~cbind(.x, .y[,1:(ncol(.y)-1)]))
map(modSep_prop, head)

## save
save(modSep, modSep_ilr, file = "differential_models_ilr.Rdata")

#### Gene ####
message("\nGENE")
rse_gene_sep <- map(regions, ~rse_gene[,rse_gene$BrainRegion == .x])
map(rse_gene_sep, dim)
# $sacc
# [1] 25212   551
# 
# $amyg
# [1] 25212   540

outGene <- map2(rse_gene_sep, modSep, ~run_DE(rse = .x, model = .y, save_eBayes = TRUE))
outGene_ilr <- map2(rse_gene_sep, modSep_ilr, ~run_DE(rse = .x, model = .y, save_eBayes = TRUE))
outGene_prop <- map2(rse_gene_sep, modSep_prop, ~run_DE(rse = .x, model = .y, save_eBayes = TRUE))

## Extract vals
ebGene <- list(outGene, outGene_ilr, outGene_prop)
ebGene <- map(ebGene, ~map(.x, "eBayes"))
names(ebGene) <- list("no_deconvo","ilr","prop")

ebGene <- transpose(ebGene)
names(ebGene)

## Extract  outGene
outGene <- map(outGene, "topTable")
outGene_ilr <- map(outGene_ilr, "topTable")
outGene_prop <- map(outGene_prop, "topTable")

## get topTable stats for 1 coef at a time
dx_coef <- list(ctrl = "PrimaryDxControl", bp = "PrimaryDxBipolar")

outGene_single_coef <- map(ebGene, function(region){
  tt_region <- map(dx_coef, function(dx){
        tt_model <- map(region, ~topTable(.x, coef= dx, p.value = 1, number=nrow(rse_gene), sort.by = "none"))
    return(tt_model)
    })
  return(tt_region)
}
)

save(outGene, outGene_ilr, outGene_single_coef,file = here("differential_expression","data", "qSVA_MDD_gene_DEresults_ilr.rda"))

#### Exon ####
message("\nEXON")
load(here('exprs_cutoff','rse_exon.Rdata'), verbose=TRUE)
rse_exon <- map(regions, ~rse_exon[,rse_exon$BrainRegion == .x])
map(rse_exon, dim)

outExon <- map2(rse_exon, modSep, run_DE)
save(outExon, file = here("differential_expression","data", "qSVA_MDD_exon_DEresults.rda"))
rm(rse_exon)

#### Jxn ####
message("\nJXN")
load(here('exprs_cutoff','rse_jxn.Rdata'), verbose=TRUE)
rse_jxn <- map(regions, ~rse_jxn[,rse_jxn$BrainRegion == .x])
map(rse_jxn, dim)

outjxn <- map2(rse_jxn, modSep, run_DE)
save(outjxn, file = here("differential_expression","data", "qSVA_MDD_jxn_DEresults.rda"))
rm(rse_jxn)

#### Tx ####
message("\nTX")
load(here('exprs_cutoff','rse_tx.Rdata'), verbose=TRUE)
rse_tx <- map(regions, ~rse_tx[,rse_tx$BrainRegion == .x])
map(rse_tx, dim)

outtx <- map2(rse_tx, modSep, ~run_DE(.x, .y, run_voom = FALSE))
save(outtx, file = here("differential_expression","data", "qSVA_MDD_tx_DEresults.rda"))


#sgejobs::job_single('qSV_model_DE_analysis_ilr', create_shell = TRUE, memory = '80G', command = "Rscript qSV_model_DE_analysis_ilr.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
