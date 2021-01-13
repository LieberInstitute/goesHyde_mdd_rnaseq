
library(jaffelab)
library(SummarizedExperiment)
library(variancePartition)
library(compositions)
library(purrr)
library(sessioninfo)
library(here)

#load objects
load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)
pd <- as.data.frame(colData(rse_gene))
# pd$Sex <- pd$Sex == "M"
 
modJoint = ~(1|PrimaryDx)*(1|BrainRegion) + AgeDeath + (1|Sex) + 
                          snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
                          mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr

#### Var Partion with Model terms ####
# form <- ~AgeDeath + (1|BrainRegion) + (1|Sex) + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5+ 
#   mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr

varPart <- fitExtractVarPartModel(exprObj =assays(rse_gene)$counts, formula = modJoint, data = pd)
vp <- sortCols(varPart)
vp_violin <- plotVarPart( vp )

ggsave(plot = vp_violin, filename = "plots/vp_violin.png")

#### Var Partition with cell type props ####
load(here("deconvolution","data","est_prop_top5.Rdata"), verbose = TRUE)

est_prop_ilr <- map(est_prop, ~ilr(.x$Est.prop.weighted))
est_prop_ilr <- map(est_prop_ilr, function(x){
  colnames(x) <- paste0("cfV",1:ncol(x))
  return(x)
})

## sACC
rse_gene_sacc <- rse_gene[,rse_gene$BrainRegion == "sACC"]
pd_sacc <- cbind(as.data.frame(colData(rse_gene_sacc)),est_prop_ilr$sacc_specific)

mod_sacc = ~(1|PrimaryDx) + AgeDeath + (1|Sex) + 
  snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
  mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr +
  cfV1 + cfV2 + cfV3 + cfV4 + cfV5 + cfV6 + cfV7 + cfV8 +cfV9 

message("VarPart sACC")
varPart_sacc <- fitExtractVarPartModel(exprObj =assays(rse_gene_sacc)$counts, formula = mod_sacc, data = pd_sacc)
vp_sacc <- sortCols(varPart_sacc)
vp_violin <- plotVarPart(vp_sacc)

ggsave(plot = vp_violin, filename = "plots/vp_violin_sacc.png")

## Amygdala
rse_gene_amyg <- rse_gene[,rse_gene$BrainRegion == "Amygdala"]
pd_amyg <- cbind(as.data.frame(colData(rse_gene_amyg)),est_prop_ilr$amyg_specific)

mod_amyg = ~(1|PrimaryDx) + AgeDeath + (1|Sex) + 
  snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
  mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr +
  cfV1 + cfV2 + cfV3 + cfV4 + cfV5 + cfV6 + cfV7 + cfV8

message("VarPart Amyg")
varPart_amyg <- fitExtractVarPartModel(exprObj =assays(rse_gene_amyg)$counts, formula = mod_amyg, data = pd_amyg)
vp_amyg <- sortCols(varPart_amyg)
vp_violin <- plotVarPart(vp_amyg)

ggsave(plot = vp_violin, filename = "plots/vp_violin_amyg.png")

save(varPart, varPart_sacc, varPart_amyg, file = "varPart.Rdata")

# sgejobs::job_single('variancePartition_exploration', create_shell = TRUE, queue= 'bluejay', memory = '150G', command = "Rscript variancePartition_exploration.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
