library(RColorBrewer)
library(ggrepel)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(jaffelab)
library(here)
library(MuSiC)
library(Biobase)
library(xbioc)
library(dplyr)
library(reshape2)

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
sce.hpc$cellType.Broad <- ss(as.character(sce.hpc$cellType.split), "\\.", 1)
table(sce.hpc$cellType.Broad)
# Astro Excit Inhib Micro Oligo   OPC Tcell 
# 1343   602   370  1253  5885   864    26


## Load rse_gene data
load(here("data", "rse_gene_GoesZandi.rda"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$Symbol
## Create input Expression sets for MuSiC
rse_gene_sACC <- rse_gene[!duplicated(rownames(rse_gene)),rse_gene$BrainRegion == "sACC"]
es_gene_sACC <- ExpressionSet(assayData = assays(rse_gene_sACC)$counts)

#create unique colnames
sce.hpc$uniqueID <- paste0(sce.hpc$donor, "_", sce.hpc$Barcode)
colnames(sce.hpc) <- sce.hpc$uniqueID
rownames(sce.hpc) <- rowData(sce.hpc)$Symbol
sce.hpc <- sce.hpc[!duplicated(rownames(sce.hpc)),]
# create pheno Data
pd_sce_sacc <- as.data.frame(colData(sce.hpc)[c("cellType.Broad", "uniqueID","sizeFactor","sizeFactor")])
# create 
es_sc_sACC <- ExpressionSet(assayData = as.matrix(assays(sce.hpc)$counts),
                            phenoData=AnnotatedDataFrame(pd_sce_sacc))


table(rownames(es_sc_sACC) %in% rownames(es_gene_sACC))
# FALSE  TRUE 
# 11120 22394 

est_prop_sacc = music_prop(bulk.eset = es_gene_sACC, 
                           sc.eset = es_sc_sACC, 
                           clusters = 'cellType.Broad',
                           samples = 'uniqueID')
save(est_prop_sacc, file = "prop_sacc.Rdata")

head(est_prop_sacc$Est.prop.weighted)
# Oligo     Micro         OPC       Inhib     Astro       Excit     Tcell
# R14149_psychENCODE_MDD 0.1798537 0.1729239 0.000000000 0.028646404 0.2713655 0.079067756 0.2681427
# R14152_psychENCODE_MDD 0.1227360 0.1772932 0.000976527 0.027257837 0.3109532 0.077233664 0.2835495
# R14153_psychENCODE_MDD 0.1002896 0.1692204 0.000000000 0.037222169 0.2556996 0.077744479 0.3598238
# R14178_psychENCODE_MDD 0.1726896 0.1882805 0.000000000 0.030266292 0.2606221 0.067435932 0.2807056
# R14179_psychENCODE_MDD 0.5888944 0.0694142 0.000000000 0.007238962 0.1807860 0.003551994 0.1501144
# R14218_psychENCODE_MDD 0.1785658 0.1046717 0.000000000 0.025534163 0.3298356 0.056813870 0.3045789
load("prop_sacc.Rdata", verbose = TRUE)

cell_bias <- music_basis(es_sc_sACC,
                   clusters = 'cellType.Broad',
                   samples = 'uniqueID')
# cell_bias$M.S
# Oligo     Micro       OPC     Inhib     Astro     Excit     Tcell 
# 6770.924  4980.110  9747.517 33398.770  7848.042 33804.591  3095.385 

pd_sacc <- as.data.frame(colData(rse_gene_sACC))
pd_sacc <- cbind(pd_sacc, est_prop_sacc$Est.prop.weighted)

## plot 
bp_dx <- pd_sacc %>% select(PrimaryDx, names(cell_bias$M.S)) %>%
  melt(id.vars = "PrimaryDx") %>%
  rename(`Cell Type` = variable, Prop = value) %>%
  ggplot(aes(`Cell Type`, Prop, fill = `PrimaryDx`)) +
  geom_boxplot() +
  labs(title = "Distribution of Cell Type Propotions",
       subtitle = "sACC samples - MuSiC defaults")

ggsave(filename = "sACC_CellType_boxplot.png", plot = bp_dx, width = 10)
