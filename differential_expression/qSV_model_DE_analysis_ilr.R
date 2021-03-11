#####################################
# from /dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/R_scripts/DEG_scripts/qSV_model_analysis_DLPFC_PTSD.R
#####################################

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)
library(edgeR)
library(tidyverse)
library(sessioninfo)
library(here)

source("run_DE.R")

##### Load Data ####

## load rse
load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)
pd = colData(rse_gene)

## qSV data
load(here("differential_expression","data","qSV_mat.Rdata"), verbose = TRUE)

## load deconvolution data
load(here("deconvolution","data","est_prop_MuSiC.Rdata"), verbose = TRUE)
load(here("deconvolution","data","est_prop_Bisque.Rdata"), verbose = TRUE)
est_prop <- list(music =  map(est_prop_music[c("sacc_specific","amyg_specific")], "Est.prop.weighted"),
                 bisque = map(est_prop_bisque[c("sacc_specific","amyg_specific")], "bulk.props"))

ilr <- list(music =  map(est_prop_music[c("sacc_specific","amyg_specific")], "ilr"),
                 bisque = map(est_prop_bisque[c("sacc_specific","amyg_specific")], "ilr"))

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
modSep_ilr <- map(ilr, function(ilr_method){
  map2(modSep, ilr_method, ~cbind(.x, .y))
})
map(modSep_ilr, ~map(.x,colnames))

## Add prop terms
map(est_prop, colnames)

modSep_prop <- map(est_prop, function(prop_method){
  map2(modSep, prop_method, ~cbind(.x, .y[,1:ncol(.y)-1]))
})
map(modSep_prop, ~map(.x,colnames))

## save
save(modSep, modSep_ilr, modSep_prop, file = "differential_models_ilr.Rdata")

#### Gene ####
message("\nGENE")
rse_gene_sep <- map(regions, ~rse_gene[,rse_gene$BrainRegion == .x])
map(rse_gene_sep, dim)
# $sacc
# [1] 25212   551
# 
# $amyg
# [1] 25212   540

## no deconvlution model
outGene <- map2(rse_gene_sep, modSep, ~run_DE(rse = .x, model = .y, save_eBayes = TRUE))
## ilr models
outGene_ilr_music <- map2(rse_gene_sep, modSep_ilr$music, ~run_DE(rse = .x, model = .y, save_eBayes = TRUE))
outGene_ilr_bisque <- map2(rse_gene_sep, modSep_ilr$bisque, ~run_DE(rse = .x, model = .y, save_eBayes = TRUE))

## prop models
outGene_prop_music <- map2(rse_gene_sep, modSep_prop$music, ~run_DE(rse = .x, model = .y, save_eBayes = TRUE))
outGene_prop_bisque <- map2(rse_gene_sep, modSep_prop$bisque, ~run_DE(rse = .x, model = .y, save_eBayes = TRUE))

save(outGene, outGene_ilr_music, outGene_ilr_bisque, outGene_prop_music, outGene_prop_bisque, 
     file = here("differential_expression","data","outGene.Rdata"))

# load(here("differential_expression","data","outGene.Rdata"), verbose = TRUE)

## Extract vals
ebGene <- list(outGene, outGene_ilr_music, outGene_ilr_bisque, outGene_prop_music, outGene_prop_bisque)
names(ebGene) <- list("no_deconvo","ilr_music", "ilr_bisque","prop_music","prop_bisque")

gene_topTables <- map(ebGene, ~map(.x, "topTable"))
gene_topTables <- transpose(gene_topTables)

ebGene <- map(ebGene, ~map(.x, "eBayes"))
ebGene <- transpose(ebGene)

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

save(gene_topTables, outGene_single_coef, file = here("differential_expression","data","qSVA_MDD_gene_DEresults.rda"))
#### get most significant deconvo terms ####
get_top_terms <- function(eBayes, terms){
  
  tt_terms <- map_dfr(terms, function(t) {
    tt <- topTable(eBayes, coef= t, p.value = 1, number=nrow(rse_gene), sort.by = "none")
    tt$coef <- t
    return(tt)
  })
  
  tt_terms_sig <- tt_terms %>% group_by(gencodeID) %>%
    arrange(adj.P.Val) %>%
    slice(1) %>%
    ungroup()
  
  return(tt_terms_sig)
}

## define terms
prop_terms <- map(est_prop, function(ep){map(ep, ~colnames(.x)[1:ncol(.x) - 1])})
prop_terms <- transpose(prop_terms)

prop_coef_explore <- map2(ebGene, prop_terms, function(eb_region, terms_region){
  eb_region <- eb_region[grepl("prop", names(eb_region))]
  map2(eb_region, terms_region, ~get_top_terms(.x, .y)%>%
         rename(coef_prop = coef))
})

ilr_terms <- map(ilr, ~map(.x, colnames))
ilr_terms <- transpose(ilr_terms)

ilr_coef_explore <- map2(ebGene, ilr_terms, function(eb_region, terms_region){
  eb_region <- eb_region[grepl("ilr", names(eb_region))]
  map2(eb_region, terms_region, ~get_top_terms(.x, .y) %>%
         rename(coef_ilr = coef))
})

combo_coef_explore <- map2(ilr_coef_explore, prop_coef_explore, function(ilr_region, prop_region){
   combo_method <- map2(ilr_region, prop_region, ~.x %>% 
    select("gencodeID", "coef_ilr") %>% 
    left_join(.y %>% select("gencodeID", "coef_prop")))
  # message(names(combo_method))
  combo_method2 <- map2(combo_method, names(combo_method), ~add_column(.x, "method" = ss(.y,"_",2)))
  combo_region <- do.call("rbind", combo_method2)
  return(combo_region)
} )

map(combo_coef_explore, ~count(.x, method, coef_ilr, coef_prop) %>% arrange(-n))
save(prop_coef_explore, ilr_coef_explore, combo_coef_explore, file = here("differential_expression","data", "deconvo_coef_explore.Rdata"))

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
