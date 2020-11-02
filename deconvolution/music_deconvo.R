
library(SummarizedExperiment)
library(SingleCellExperiment)
library(jaffelab)
library(here)
library(MuSiC)
library(Biobase)
library(xbioc)

## Load rse_gene data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

#### sacc Data ####
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose = TRUE)
dim(sce.sacc)
# [1] 33538  7047
## Exclude Ambig.lowNtrxts
sce.sacc <- sce.sacc[,sce.sacc$cellType != "Ambig.lowNtrxts",]

## Add cellType.broad
sce.sacc$cellType.Broad <- ss(as.character(sce.sacc$cellType), "\\.", 1)
table(sce.sacc$cellType.Broad)
# Ambig Astro Excit Inhib Micro Oligo   OPC 
# 43   623  1266   830   502  3252   531 

## Match rownames
rownames(sce.sacc) <- rowData(sce.sacc)$ID
message("Matching genes:")
table(rownames(sce.sacc) %in% rownames(rse_gene))
# FALSE  TRUE 
# 15021 18517 

## Create input Expression sets for MuSiC
rse_gene_sacc <- rse_gene[,rse_gene$BrainRegion == "sACC"]
es_gene_sacc <- ExpressionSet(assayData = assays(rse_gene_sacc)$counts)

#create unique colnames
sce.sacc$uniqueID <- paste0(sce.sacc$donor, "_", sce.sacc$Barcode)
colnames(sce.sacc) <- sce.sacc$uniqueID

# create pheno Data
pd_sce_sacc <- as.data.frame(colData(sce.sacc)[c("cellType","cellType.Broad", "uniqueID")])
# create 
es_sc_sacc <- ExpressionSet(assayData = as.matrix(assays(sce.sacc)$counts),
                            phenoData=AnnotatedDataFrame(pd_sce_sacc))


#### Amyg Dat ####
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_Amyg-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose = TRUE)
dim(sce.amy)
# [1] 33538  6632
table(sce.amy$cellType.split)

## Exclude Ambig.lowNtrxts
sce.amy <- sce.amy[,sce.amy$cellType.split != "Ambig.lowNtrxts",]

## Add cellType.broad
sce.amy$cellType.Broad <- ss(as.character(sce.amy$cellType.split), "\\.", 1)
table(sce.amy$cellType.Broad)
# Astro Excit Inhib Micro Oligo   OPC 
# 852   429   437   764  3473   627 

## Match rownames
rownames(sce.amy) <- rowData(sce.amy)$ID
message("Matching genes:")
table(rownames(sce.amy) %in% rownames(rse_gene))
# FALSE  TRUE 
# 15021 18517 

## Create input Expression sets for MuSiC
rse_gene_amyg <- rse_gene[,rse_gene$BrainRegion == "Amygdala"]
es_gene_amyg <- ExpressionSet(assayData = assays(rse_gene_amyg)$counts)

#create unique colnames
sce.amy$uniqueID <- paste0(sce.amy$donor, "_", sce.amy$Barcode)
colnames(sce.amy) <- sce.amy$uniqueID

# create pheno Data
pd_sce_amyg <- as.data.frame(colData(sce.amy)[c("cellType.split","cellType.Broad", "uniqueID")])
# create 
es_sc_amyg <- ExpressionSet(assayData = as.matrix(assays(sce.amy)$counts),
                            phenoData=AnnotatedDataFrame(pd_sce_amyg))


#### estimate cell type props ####

message("**** MUSIC for sACC ****")
est_prop_sacc = music_prop(bulk.eset = es_gene_sacc, 
                           sc.eset = es_sc_sacc, 
                           clusters = 'cellType',
                           samples = 'uniqueID')
save(est_prop_sacc, file = "prop_sacc.Rdata")

est_prop_broad_sacc = music_prop(bulk.eset = es_gene_sacc, 
                                 sc.eset = es_sc_sacc, 
                                 clusters = 'cellType.Broad',
                                 samples = 'uniqueID')
save(est_prop_broad_sacc, file = "prop_broad_sacc.Rdata")

message("**** MUSIC for Amyg ****")
est_prop_amyg = music_prop(bulk.eset = es_gene_amyg, 
                           sc.eset = es_sc_amyg, 
                           clusters = 'cellType.split',
                           samples = 'uniqueID')
save(est_prop_amyg, file = "prop_amyg.Rdata")

est_prop_broad_amyg = music_prop(bulk.eset = es_gene_amyg, 
                                 sc.eset = es_sc_amyg, 
                                 clusters = 'cellType.Broad',
                                 samples = 'uniqueID')
save(est_prop_broad_amyg, file = "prop_broad_amyg.Rdata")




#### Top 40 data ####
top40_amyg <- read.csv("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/tables/top40genesLists_Amyg-n2_cellType.split_SN-LEVEL-tests_May2020.csv")
top40_sacc <- read.csv("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/tables/top40genesLists_sacc-n2_cellType_SN-LEVEL-tests_May2020.csv")

## sacc
top40all_sacc <- unique(unlist(top40_sacc[,grepl("1vAll", colnames(top40_sacc))]))
length(top40all_sacc)
# [1] 388
top40all_ensm_sacc <-  rowData(rse_gene)$ensemblID[rowData(rse_gene)$Symbol %in% top40all_sacc]
top40all_ensm_sacc <- top40all_ensm_sacc[top40all_ensm_sacc %in% rownames(es_sc_sacc)]
length(top40all_ensm_sacc)
# [1] 332
## estimate cell type props
est_prop_top40_sacc = music_prop(bulk.eset = es_gene_sacc[top40all_ensm_sacc,], 
                                 sc.eset = es_sc_sacc[top40all_ensm_sacc,], 
                                 clusters = 'cellType',
                                 samples = 'uniqueID')
save(est_prop_top40_sacc, file = "prop_top40_sacc.Rdata")

##Amyg
top40all_amyg <- unique(unlist(top40_amyg[,grepl("1vAll", colnames(top40_amyg))]))
length(top40all_amyg)
# [1] 441
top40all_ensm_amyg <-  rowData(rse_gene)$ensemblID[rowData(rse_gene)$Symbol %in% top40all_amyg]
top40all_ensm_amyg <- top40all_ensm_amyg[top40all_ensm_amyg %in% rownames(es_sc_amyg)]
## estimate cell type props
est_prop_top40_amyg = music_prop(bulk.eset = es_gene_sacc[top40all_ensm_amyg,], 
                                 sc.eset = es_sc_sacc[top40all_ensm_amyg,], 
                                 clusters = 'cellType',
                                 samples = 'uniqueID')
save(est_prop_top40_amyg, file = "prop_top40_amyg.Rdata")

# sgejobs::job_single('music_deconvo', create_shell = TRUE, queue= 'bluejay', memory = '100G', command = "Rscript music_deconvo.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


