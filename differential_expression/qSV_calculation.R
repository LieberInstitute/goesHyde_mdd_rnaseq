#####################################
# from /dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/R_scripts/DEG_scripts/qSV_model_analysis_DLPFC_PTSD.R
#####################################

library(jaffelab)
library(SummarizedExperiment)
library(sva)
# library(readxl)
# library(devtools)
# library(edgeR)
library(sessioninfo)
library(here)

##### Load rse data, examine ####

#load objects
load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)
pd = colData(rse_gene)

table(pd$Experiment)
# psychENCODE_MDD  psychENCODE_BP 
# 588             503 

## load degradation data
load(here("data","degradation_rse_MDDseq_BiPSeq_BothRegions.Rdata"), verbose = TRUE)

##### get qSVs ####
modJoint = model.matrix(~PrimaryDx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 +
	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr,
	data=colData(rse_gene))

colnames(modJoint)

#  [1] "(Intercept)"                      "PrimaryDxControl"
#  [3] "PrimaryDxBipolar"                 "BrainRegionsACC"
#  [5] "AgeDeath"                         "SexM"
#  [7] "snpPC1"                           "snpPC2"
#  [9] "snpPC3"                           "mitoRate"
# [11] "rRNA_rate"                        "totalAssignedGene"
# [13] "RIN"                              "PrimaryDxControl:BrainRegionsACC"
# [15] "PrimaryDxBipolar:BrainRegionsACC"

degExprs = log2(assays(cov_rse)$count+1)
k = num.sv(degExprs, modJoint)
message("k=", k)
# k=26
qSV_mat = prcomp(t(degExprs))$x[,1:k]
save(qSV_mat, file = "qSV_mat.Rdata")
varExplQsva = jaffelab::getPcaVars(prcomp(t(degExprs)))
varExplQsva[1:k]

# [1] 67.700  4.860  3.010  1.990  1.430  1.210  1.120  0.803  0.778  0.677
# [11]  0.574  0.521  0.487  0.444  0.347  0.342  0.329  0.300  0.246  0.237
# [21]  0.230  0.223  0.205  0.196  0.190  0.174


sum(varExplQsva[1:k]) # 88.623%
# [1] 88.125


#sgejobs::job_single('qSV_calculation', create_shell = TRUE, memory = '80G', command = "Rscript qSV_calculations.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
