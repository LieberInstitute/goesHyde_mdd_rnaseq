
library(jaffelab)
library(SummarizedExperiment)
library(variancePartition)
library(sessioninfo)
library(here)

#load objects
load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)
pd <- colData(rse_gene)
pd$Sex <- pd$Sex == "M"
  
load(here("deconvolution","data","est_prop_top5.Rdata"))

# modJoint = model.matrix(~PrimaryDx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 +
#                           mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr,
#                         data=colData(rse_gene))

form <- ~AgeDeath + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5+ 
  mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr


varPart <- fitExtractVarPartModel(exprObj =assays(rse_gene)$counts, formula = form, data = pd)

