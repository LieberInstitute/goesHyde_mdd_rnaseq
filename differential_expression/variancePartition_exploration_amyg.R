
library(jaffelab)
library(SummarizedExperiment)
library(variancePartition)
library(compositions)
library(purrr)
library(limma)
library(edgeR)
library(sessioninfo)
library(here)

# load("varPart.Rdata", verbose = TRUE)

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

## Amygdala
rse_gene_amyg <- rse_gene[,rse_gene$BrainRegion == "Amygdala"]
pd_amyg <- cbind(as.data.frame(colData(rse_gene_amyg)),est_prop_ilr$amyg_specific)

mod_amyg <- model.matrix(~PrimaryDx, pd_amyg)
gExpr_amyg <- calcNormFactors(rse_gene_amyg)
vobjGenes_amyg <- voom(gExpr_amyg, mod_amyg)

mod_amyg = ~(1|PrimaryDx) + AgeDeath + (1|Sex) + 
  snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
  mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr +
  ilr_1 + ilr_2 + ilr_3 + ilr_4 + ilr_5 + ilr_6 + ilr_7 + ilr_8

message("VarPart Amyg")
varPart_amyg <- fitExtractVarPartModel(exprObj = vobjGenes_amyg, formula = mod_amyg, data = pd_amyg)
vp_amyg <- sortCols(varPart_amyg)
vp_violin <- plotVarPart(vp_amyg) + labs(title = "mod_amyg")

ggsave(plot = vp_violin, filename = "plots/vp_violin_amyg.png", width = 10)

save(varPart_amyg, file = "varPart_amyg.Rdata")

# sgejobs::job_single('variancePartition_exploration_amyg', create_shell = TRUE, queue= 'bluejay', memory = '150G', command = "Rscript variancePartition_exploration_amyg.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
