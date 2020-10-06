
library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(edgeR)
library(limma)
library(sessioninfo)
library(here)

load(here("data","rse_gene_overlap_GoesZandi.rda"), verbose = TRUE)

mod <- model.matrix(~Experiment + BrNum, data = colData(rse_gene_ol))
colnames(model)

dge = DGEList(counts = assays(rse_gene_ol)$counts, genes = rowData(rse_gene_ol))
vGene = voom(dge, mod)
### limma commads
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)

#compute T-stat
outGene_ol = topTable(eBGene,coef=2,p.value = 1,number=nrow(rse_gene_ol))

table(outGene$adj.P.Val < 0.05)
table(outGene$adj.P.Val < 0.01)

save(outGene_ol, file = "outGene_overlap.rda")

#sgejobs::job_single('experiment_DE_analysis', create_shell = TRUE, memory = '5G', command = "Rscript experiment_DE_analysis.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
