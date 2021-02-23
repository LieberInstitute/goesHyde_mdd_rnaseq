
library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(compositions)
library(sessioninfo)
library(here)


##### Load rse data ####
#load objects
load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)
pd = as.data.frame(colData(rse_gene))

#### Var Partition with cell type props ####
## load ilr data
load(here("deconvolution","data","est_prop_ilr.Rdata"), verbose = TRUE)
est_prop_broad <- rbind(est_prop_ilr$sacc_broad, est_prop_ilr$amyg_broad)
dim(est_prop_broad)
pd <- cbind(pd, est_prop_ilr)

## load degradation data
load(here("data","degradation_rse_MDDseq_BiPSeq_BothRegions.Rdata"), verbose = TRUE)

pd <- cbind(pd,est_prop_ilr)
##### get qSVs ####
modJoint = model.matrix(~PrimaryDx*BrainRegion + AgeDeath + Sex + 
                          snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
                          mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr +
                          ilr_1 + ilr_2 + ilr_3 + ilr_4 + ilr_5,
                        data=pd)
colnames(modJoint)

# [1] "(Intercept)"                      "PrimaryDxControl"                 "PrimaryDxBipolar"                
# [4] "BrainRegionsACC"                  "AgeDeath"                         "SexM"                            
# [7] "snpPC1"                           "snpPC2"                           "snpPC3"                          
# [10] "snpPC4"                           "snpPC5"                           "mitoRate"                        
# [13] "rRNA_rate"                        "totalAssignedGene"                "RIN"                             
# [16] "ERCCsumLogErr"                    "ilr_1"                            "ilr_2"                           
# [19] "ilr_3"                            "ilr_4"                            "ilr_5"                           
# [22] "PrimaryDxControl:BrainRegionsACC" "PrimaryDxBipolar:BrainRegionsACC"

degExprs = log2(assays(cov_rse)$count+1)
k = num.sv(degExprs, modJoint)
message("k=", k)
# k=26 - model with no ilr
# k=25
# qSV_mat = prcomp(t(degExprs))$x[,1:k]
# 
# ## Save qSVmat
# save(qSV_mat, file = "qSV_mat.Rdata")

## Explore variance explained
varExplQsva = jaffelab::getPcaVars(prcomp(t(degExprs)))
varExplQsva[1:k]

# [1] 67.700  4.860  3.010  1.990  1.430  1.210  1.120  0.803  0.778  0.677  0.574  0.521  0.487  0.444  0.347  0.342
# [17]  0.329  0.300  0.246  0.237  0.230  0.223  0.205  0.196  0.190

sum(varExplQsva[1:k])
# [1] 88.449

#sgejobs::job_single('qSV_exploration', create_shell = TRUE, memory = '80G', command = "Rscript qSV_exploration.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
