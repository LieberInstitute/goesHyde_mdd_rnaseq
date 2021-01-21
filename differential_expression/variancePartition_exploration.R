
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
pd <- as.data.frame(colData(rse_gene))

mod = model.matrix(~PrimaryDx*BrainRegion, pd)
gExpr <- calcNormFactors(rse_gene)
vobjGenes <- voom(gExpr, mod)

#### Var Partition with cell type props ####
load(here("deconvolution","data","est_prop_top5.Rdata"), verbose = TRUE)

## compute ilr for broad cell types
est_prop_broad <- rbind(est_prop$sacc_broad$Est.prop.weighted, est_prop$amyg_broad$Est.prop.weighted)
dim(est_prop_broad)
est_prop_ilr <- ilr(est_prop_broad)
colnames(est_prop_ilr) <-paste0("ilr_",1:ncol(est_prop_ilr))

pd <- cbind(pd, est_prop_ilr)
## Linear 
formJoint <- list( linear = ~PrimaryDx*BrainRegion + AgeDeath + Sex +
                     snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
                     mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr,
                   
                   mixed = ~(1|PrimaryDx) + (1|BrainRegion) + (1|PrimaryDx:BrainRegion) + AgeDeath + (1|Sex) +
                     snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
                     mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr,
                   
                   linear_ilr = ~PrimaryDx*BrainRegion + AgeDeath + Sex +
                     snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
                     mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr+
                     ilr_1 + ilr_2 + ilr_3 + ilr_4 + ilr_5,
                   
                   mixed_ilr = ~(1|PrimaryDx) + (1|BrainRegion) + (1|PrimaryDx:BrainRegion) + AgeDeath + (1|Sex) +
                     snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
                     mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr +
                     ilr_1 + ilr_2 + ilr_3 + ilr_4 + ilr_5)


#### Var Partion with linear Model ####
varPart <- map(formJoint, ~fitExtractVarPartModel(exprObj =vobjGenes, formula = .x, data = pd))
vp <- map(varPart, sortCols)

vp_violin <- map2(vp, names(vp), ~plotVarPart(.x) + labs(title = paste("modJoint -",.y)))

walk2(vp_violin, names(vp_violin), ggsave(plot = vp_violin, filename = paste0("plots/vp_violin_",.y,".png")))
save(varPart, formJoint, file = "varPart.Rdata")


# sgejobs::job_single('variancePartition_exploration', create_shell = TRUE, queue= 'bluejay', memory = '150G', command = "Rscript variancePartition_exploration.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
