#####################################
# from /dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/R_scripts/DEG_scripts/qSV_model_analysis_DLPFC_PTSD.R
#####################################

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)
library(edgeR)
library(here)
library(RColorBrewer)

#load objects
load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)
load(here('exprs_cutoff','rse_exon.Rdata'), verbose=TRUE)
load(here('exprs_cutoff','rse_jxn.Rdata'), verbose=TRUE)
load(here('exprs_cutoff','rse_tx.Rdata'), verbose=TRUE)


## load degradation data
load(here("data","degradation_rse_MDDseq_BiPSeq_BothRegions.Rdata"), verbose = TRUE)

## load models
load("differental_models.Rdata", verbose = TRUE)

##### get qSVs ####
degExprs = log2(assays(cov_rse)$count+1)
k = num.sv(degExprs, modJoint) # 26
qSV_mat = prcomp(t(degExprs))$x[,1:k]
varExplQsva = getPcaVars(prcomp(t(degExprs)))
varExplQsva[1:k]
sum(varExplQsva[1:k]) # 88.623%

## load DE results
load("qSVA_MDD_gene_DEresults.rda", verbose=TRUE)
load("qSVA_MDD_exon_DEresults.rda", verbose=TRUE)
load("qSVA_MDD_jxn_DEresults.rda", verbose=TRUE)
load("qSVA_MDD_tx_DEresults.rda", verbose=TRUE)

load("outGene_overlap.rda", verbose = TRUE)

#### expression ####
rse_gene$Dx = as.factor(rse_gene$PrimaryDx)
rse_gene$Dx = factor(rse_gene$Dx, levels(rse_gene$Dx)[c(2,3,1)] )

pdSacc = colData(rse_gene)[sACC_Index,]
pdAmyg = colData(rse_gene)[Amyg_Index,]

gRpkm = recount::getRPKM(rse_gene,"Length")
eRpkm = recount::getRPKM(rse_exon,"Length")
jRpkm = recount::getRPKM(rse_jxn,"Length")
tRpkm = assays(rse_tx)$tpm

gExprs = as.matrix(log2(gRpkm+1))
eExprs = as.matrix(log2(eRpkm+1))
jExprs = as.matrix(log2(jRpkm+1))
tExprs = as.matrix(log2(tRpkm+1))

gSaccExprs = cleaningY(gExprs[,sACC_Index], mod_sACC, P=3)
gAmygExprs = cleaningY(gExprs[,Amyg_Index], mod_Amyg, P=3)

eSaccExprs = cleaningY(eExprs[,sACC_Index], mod_sACC, P=3)
eAmygExprs = cleaningY(eExprs[,Amyg_Index], mod_Amyg, P=3)

jSaccExprs = cleaningY(jExprs[,sACC_Index], mod_sACC, P=3)
jAmygExprs = cleaningY(jExprs[,Amyg_Index], mod_Amyg, P=3)

tSaccExprs = cleaningY(tExprs[,sACC_Index], mod_sACC, P=3)
tAmygExprs = cleaningY(tExprs[,Amyg_Index], mod_Amyg, P=3)


#### plot residualize expression ####
source("residExp_beeswarm.R")

## sACC
topGenes_residExp_beeswarm(outGene_sACC, gSaccExprs, pdSacc, "plots/top_genes_MDD_vs_cnt_sACC_beeSwarm.pdf")
topGenes_residExp_beeswarm(outExon_sACC, eSaccExprs, pdSacc, "plots/top_exons_MDD_vs_cnt_sACC_beeSwarm.pdf")
topGenes_residExp_beeswarm(outJxn_sACC, jSaccExprs, pdSacc, "plots/top_jxns_MDD_vs_cnt_sACC_beeSwarm.pdf")
topGenes_residExp_beeswarm(outTx_sACC, tSaccExprs, pdSacc, "plots/top_tx_MDD_vs_cnt_sACC_beeSwarm.pdf")

## Amyg
topGenes_residExp_beeswarm(outGene_Amyg, gAmygExprs, pdAmyg, "plots/top_genes_MDD_vs_cnt_Amyg_beeSwarm.pdf")
topGenes_residExp_beeswarm(outExon_Amyg, eAmygExprs, pdAmyg, "plots/top_exons_MDD_vs_cnt_Amyg_beeSwarm.pdf")
topGenes_residExp_beeswarm(outJxn_Amyg, jAmygExprs, pdAmyg, "plots/top_jxns_MDD_vs_cnt_Amyg_beeSwarm.pdf")
topGenes_residExp_beeswarm(outTx_Amyg, tAmygExprs, pdAmyg, "plots/top_tx_MDD_vs_cnt_Amyg_beeSwarm.pdf")

#sgejobs::job_single('qSV_model_DE_plot', create_shell = TRUE, memory = '80G', command = "Rscript qSV_model_DE_plot.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
