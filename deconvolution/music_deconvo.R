library(RColorBrewer)
library(ggrepel)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(jaffelab)

load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_HPC-n3_cleaned-combined_SCE_MNTFeb2020.rda", verbose = TRUE)
dim(sce.hpc)
# [1] 33538 10444

## drop the "Ambig.lowNtrxts" from sce.hpc$cellType.split
table(sce.hpc$cellType.split ==  "Ambig.lowNtrxts")
# FALSE  TRUE 
# 10343   101
sce.hpc <- sce.hpc[,sce.hpc$cellType.split !=  "Ambig.lowNtrxts"]

## Add cellType.broad
table(sce.hpc$cellType.split)
# Ambig.lowNtrxts           Astro         Excit.1         Excit.2         Excit.3         Excit.4         Excit.5 
# 0            1343             116             117             310              26              33 
# Inhib.1         Inhib.2         Inhib.3         Inhib.4         Inhib.5           Micro           Oligo 
# 30              90             139              55              56            1253            5885 
# OPC           Tcell 
# 864              26 
sce.hpc$cellType.Broad <- ss(as.character(sce.hpc$cellType.split), "\\.", 1)
table(sce.hpc$cellType.Broad)
# Astro Excit Inhib Micro Oligo   OPC Tcell 
# 1343   602   370  1253  5885   864    26
