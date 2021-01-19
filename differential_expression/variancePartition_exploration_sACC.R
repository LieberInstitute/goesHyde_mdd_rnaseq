
library(jaffelab)
library(SummarizedExperiment)
library(variancePartition)
library(compositions)
library(purrr)
library(limma)
library(edgeR)
library(sessioninfo)
library(here)

#load objects
load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)

#### Var Partition with cell type props ####
load(here("deconvolution","data","est_prop_top5.Rdata"), verbose = TRUE)

## compute ilr 
est_prop_ilr <- map(est_prop, ~ilr(.x$Est.prop.weighted))
est_prop_ilr <- map(est_prop_ilr, function(x){
  colnames(x) <- paste0("ilr_",1:ncol(x))
  return(x)
})

#### sACC ####
rse_gene_sacc <- rse_gene[,rse_gene$BrainRegion == "sACC"]
pd_sacc <- cbind(as.data.frame(colData(rse_gene_sacc)),est_prop_ilr$sacc_specific)

mod_sacc <- model.matrix(~PrimaryDx, pd_sacc)
gExpr_sacc <- calcNormFactors(rse_gene_sacc)
vobjGenes_sacc <- voom(gExpr_sacc, mod_sacc)

mod_sacc = ~(1|PrimaryDx) + AgeDeath + (1|Sex) + 
  snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
  mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr +
  ilr_1 + ilr_2 + ilr_3 + ilr_4 + ilr_5 + ilr_6 + ilr_7 + ilr_8 +ilr_9 

message("VarPart sACC")
varPart_sacc <- fitExtractVarPartModel(exprObj = vobjGenes_sacc, formula = mod_sacc, data = pd_sacc)
vp_sacc <- sortCols(varPart_sacc)
vp_violin <- plotVarPart(vp_sacc)+ labs(title = "mod_sacc")

ggsave(plot = vp_violin, filename = "plots/vp_violin_sACC.png", width = 10)
save(varPart_sacc, file = "varPart_sacc.Rdata")

# sgejobs::job_single('variancePartition_exploration_sACC', create_shell = TRUE, queue= 'bluejay', memory = '150G', command = "Rscript variancePartition_exploration_sACC.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
