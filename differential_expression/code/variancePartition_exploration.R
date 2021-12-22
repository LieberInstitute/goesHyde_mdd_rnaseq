
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
load(here("differential_expression","data","qSV_mat.Rdata"), verbose = TRUE)
cbind(colData(rse_gene),qSV_mat)
#### Var Partition with cell type props ####
message("Subset data")
rse_gene <- rse_gene[,rse_gene$BrainRegion == test_region & rse_gene$PrimaryDx %in% c("Control", test_dx)]
dim(rse_gene)
pd <- as.data.frame(colData(rse_gene))

mod <- model.matrix(~PrimaryDx, pd)
gExpr <- calcNormFactors(rse_gene)
vobjGenes <- voom(gExpr, mod)

# linear mixed model where all categorical variables are modeled as random
# effects and all continuous variables are fixed effects. The function lmer
# from lme4 is used to fit this model.

form = list(Sep = ~(1|PrimaryDx) + AgeDeath + (1|Sex) + 
              snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
              mitoRate + rRNA_rate + totalAssignedGene + RIN + abs(ERCCsumLogErr) +
              PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + 
              PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + 
              PC21 + PC22 + PC23 + PC24 + PC25 + PC26,
            
            Sep_cf = ~(1|PrimaryDx) + AgeDeath + (1|Sex) + 
              snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
              mitoRate + rRNA_rate + totalAssignedGene + RIN + abs(ERCCsumLogErr) +
              Astro + Endo + Macro + Micro + Mural + Oligo + OPC + Tcell + Excit +
              PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + 
              PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + 
              PC21 + PC22 + PC23 + PC24 + PC25 + PC26)



message("VarPart Amyg")
varPart <- purrr::map(form, ~fitExtractVarPartModel(exprObj = vobjGenes, formula = .x, data = pd))
save(varPart, file = here("differential_expression", "data","variance_partition", paste0("varPart_",test_region,"_",test_dx,".Rdata")))

# sgejobs::job_single('variancePartition_exploration', create_shell = TRUE, queue= 'bluejay', memory = '150G', command = "Rscript variancePartition_exploration.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
