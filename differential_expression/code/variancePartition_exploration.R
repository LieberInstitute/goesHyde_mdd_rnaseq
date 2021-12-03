
library(jaffelab)
library(SummarizedExperiment)
library(variancePartition)
library(limma)
library(edgeR)
library(sessioninfo)
library(here)


args = commandArgs(trailingOnly=TRUE)[[1]]
test_region = ss(args,"_")
test_dx = ss(args,"_",2)

message(paste("Running", test_dx, "vs. Control -", test_region))
#load objects
load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)

#### Var Partition with cell type props ####
message("Subset data")
rse_gene <- rse_gene[,rse_gene$BrainRegion == test_region & rse_gene$PrimaryDx %in% c("Control", test_dx)]
dim(rse_gene)
pd <- as.data.frame(colData(rse_gene))

mod <- model.matrix(~PrimaryDx, pd)
gExpr <- calcNormFactors(rse_gene)
vobjGenes <- voom(gExpr, mod)

form = ~PrimaryDx + AgeDeath + Sex + 
  snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
  mitoRate + rRNA_rate + totalAssignedGene + RIN + abs(ERCCsumLogErr) 

message("VarPart Amyg")
varPart <- fitExtractVarPartModel(exprObj = vobjGenes, formula = form, data = pd)
save(varPart, file = here("differential_expression", "data","variance_partition", paste0("varPart_",test_region,"_",test_dx,".Rdata")))

# sgejobs::job_single('variancePartition_exploration', create_shell = TRUE, queue= 'bluejay', memory = '150G', command = "Rscript variancePartition_exploration.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
