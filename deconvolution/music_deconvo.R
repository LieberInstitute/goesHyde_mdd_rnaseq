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

## Load rse_gene data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

#### sACC Data ####
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose = TRUE)
dim(sce.sacc)
# [1] 33538  7047
table(sce.sacc$cellType)
# Ambig.lowNtrxts           Astro         Excit.1         Excit.2         Excit.3         Excit.4         Inhib.1 
# 43             623             586             410             185              85             515 
# Inhib.2           Micro           Oligo             OPC 
# 315             502            3252             531 
## Add cellType.broad
sce.sacc$cellType.Broad <- ss(as.character(sce.sacc$cellType), "\\.", 1)
table(sce.sacc$cellType.Broad)
# Ambig Astro Excit Inhib Micro Oligo   OPC 
# 43   623  1266   830   502  3252   531 

## Match rownames
rownames(sce.sacc) <- rowData(sce.sacc)$ID
table(rownames(sce.sacc) %in% rownames(rse_gene))
# FALSE  TRUE 
# 15021 18517 

## Create input Expression sets for MuSiC
rse_gene_sACC <- rse_gene[,rse_gene$BrainRegion == "sACC"]
es_gene_sACC <- ExpressionSet(assayData = assays(rse_gene_sACC)$counts)

#create unique colnames
sce.sacc$uniqueID <- paste0(sce.sacc$donor, "_", sce.sacc$Barcode)
colnames(sce.sacc) <- sce.sacc$uniqueID

# create pheno Data
pd_sce_sacc <- as.data.frame(colData(sce.sacc)[c("cellType.Broad", "uniqueID")])
# create 
es_sc_sACC <- ExpressionSet(assayData = as.matrix(assays(sce.sacc)$counts),
                            phenoData=AnnotatedDataFrame(pd_sce_sacc))


## estimate cell type props
est_prop_sacc = music_prop(bulk.eset = es_gene_sACC, 
                           sc.eset = es_sc_sACC, 
                           clusters = 'cellType.Broad',
                           samples = 'uniqueID')
save(est_prop_sacc, file = "prop_sacc.Rdata")

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
table(rownames(sce.amy) %in% rownames(rse_gene))
# FALSE  TRUE 
# 15021 18517 

## Create input Expression sets for MuSiC
rse_gene_amyg <- rse_gene[,rse_gene$BrainRegion == "Amyg"]
es_gene_amyg <- ExpressionSet(assayData = assays(rse_gene_amyg)$counts)

#create unique colnames
sce.amy$uniqueID <- paste0(sce.amy$donor, "_", sce.amy$Barcode)
colnames(sce.amy) <- sce.amy$uniqueID

# create pheno Data
pd_sce_sacc <- as.data.frame(colData(sce.amy)[c("cellType.Broad", "uniqueID")])
# create 
es_sc_sACC <- ExpressionSet(assayData = as.matrix(assays(sce.amy)$counts),
                            phenoData=AnnotatedDataFrame(pd_sce_sacc))


## estimate cell type props
est_prop_amyg = music_prop(bulk.eset = es_gene_sACC, 
                           sc.eset = es_sc_sACC, 
                           clusters = 'cellType.Broad',
                           samples = 'uniqueID')
save(est_prop_amyg, file = "prop_sacc.Rdata")




#### plot ####
# load("prop_sacc.Rdata", verbose = TRUE)
pd_sacc <- as.data.frame(colData(rse_gene_sACC))
pd_sacc <- cbind(pd_sacc, est_prop_sacc$Est.prop.weighted)


bp_sacc_dx <- pd_sacc %>% select(PrimaryDx, names(cell_bias$M.S)) %>%
  melt(id.vars = "PrimaryDx") %>%
  rename(`Cell Type` = variable, Prop = value) %>%
  ggplot(aes(`Cell Type`, Prop, fill = `PrimaryDx`)) +
  geom_boxplot() +
  labs(title = "Distribution of Cell Type Propotions",
       subtitle = "sACC samples - MuSiC defaults")

ggsave(filename = "sACC_CellType_boxplot.png", plot = bp_sacc_dx, width = 10)


#### check against qSV 1####
load(here("differential_expression" ,"qSV_mat.Rdata"), verbose = TRUE)
all(colnames(rse_gene)==rownames(qSV_mat))

pd_sacc <- cbind(pd_sacc, qSV_mat[rownames(pd_sacc),])

scatter_qsv1 <- pd_sacc %>% select(PrimaryDx, colnames(est_prop_sacc_noT$Est.prop.weighted), PC1) %>%
  melt(id.vars = c("PrimaryDx","PC1")) %>% 
  rename(`Cell Type` = variable, Prop = value) %>%
  ggplot(aes(PC1, Prop, color = `Cell Type`)) +
  geom_point() +
  facet_wrap(~`Cell Type`, scales = "free_y")+
  labs(title = "qSV1 vs. Cell Type Prop",
       subtitle = "sACC samples - MuSiC defaults - no T cell")

ggsave(filename = "sACC_CellType_qSV1.png", plot = scatter_qsv1, width = 14)

