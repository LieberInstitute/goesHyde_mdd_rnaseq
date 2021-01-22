
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
pd <- as.data.frame(colData(rse_gene))

mod = model.matrix(~PrimaryDx*BrainRegion, pd)
gExpr <- calcNormFactors(rse_gene)
vobjGenes <- voom(gExpr, mod)

## load ilr data
load(here("deconvolution","data","est_prop_ilr.Rdata"), verbose = TRUE)
est_prop_broad <- rbind(est_prop_ilr$sacc_broad, est_prop_ilr$amyg_broad)
dim(est_prop_broad)
pd <- cbind(pd, est_prop_broad)

#### Var Partition with cell type props ####
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
save(varPart, formJoint, file = "varPart.Rdata")
message("SAVED")
# load("varPart.Rdata", verbose = TRUE)

vp <- map(varPart, sortCols)

vp_violin <- map2(vp, names(vp), ~plotVarPart(.x) + labs(title = paste("modJoint -",.y)))

pdf("plots/variancePartition.pdf")
walk(vp_violin, print)
dev.off()

# sgejobs::job_single('variancePartition_exploration', create_shell = TRUE, queue= 'bluejay', memory = '150G', command = "Rscript variancePartition_exploration.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
