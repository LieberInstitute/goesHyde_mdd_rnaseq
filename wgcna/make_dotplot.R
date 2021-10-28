library(sva)
library(lmerTest)
library(SummarizedExperiment)
library(jaffelab)
library(WGCNA)
library(broom)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(readxl)
library(RColorBrewer)
library(sessioninfo)
library(here)

rm(list = ls())

########################
####COMBINED SAMPLE#####
########################

load("../exprs_cutoff/rse_gene.Rdata", verbose = TRUE)
load("../data/degradation_rse_MDDseq_BiPSeq_BothRegions.Rdata", verbose = TRUE)

rse_gene <- rse_gene[, rse_gene$PrimaryDx %in% c("Control", "MDD")]
rse_gene$Dx <- droplevels(rse_gene$PrimaryDx)
rse_gene$Dx<-relevel(rse_gene$Dx, "Control")

load("rdas/constructed_network_signed_bicor.rda", verbose=TRUE)
mergedColors = labels2colors(net$colors)
net$colorsLab = labels2colors(net$colors)
colorDat = data.frame(num = net$colors, col = net$colorsLab, 
	stringsAsFactors=FALSE)
colorDat$Label = paste0("ME", colorDat$num)
colorDat = colorDat[order(colorDat$num),]
colorDat = colorDat[!duplicated(colorDat$num),]
colorDat$numGenes = table(net$colors)[as.character(colorDat$num)]

dim(colorDat)

gList = split(rowData(rse_gene)$EntrezID, 
	factor(net$colorsLab,levels = colorDat$col))
gList = lapply(gList, function(x) as.character(x[!is.na(x)]))
univ = rowData(rse_gene)$EntrezID
univ = as.character(univ[!is.na(univ)])

###
go_lightcyan = enrichGO(gList$lightcyan, univ = univ,
	OrgDb = "org.Hs.eg.db", ont = "ALL",
	readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 0.05)
go_pink = enrichGO(gList$pink, univ = univ,
	OrgDb = "org.Hs.eg.db", ont = "ALL",
	readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 0.05)


pdf("dotplot/Combined_sample_WGCNA_pathway_enrichment.pdf",h=6,w=12)
dotplot(go_lightcyan, showCategory = 20, title = "Combined sample lightcyan")
dotplot(go_pink, showCategory = 20, title = "Combined sample pink")
dev.off()


########################
####sACC #####
########################
rm(list = ls())

load(here('exprs_cutoff', 'rse_gene.Rdata'), verbose = TRUE)

rse_gene <- rse_gene[, rse_gene$PrimaryDx %in% c("Control", "MDD")]
rse_gene$Dx <- droplevels(rse_gene$PrimaryDx)
rse_gene_sACC <- rse_gene[, rse_gene$BrainRegion %in% c("sACC")]
rse_gene_sACC$PrimaryDx <- droplevels(rse_gene_sACC$PrimaryDx)
rse_gene_sACC$PrimaryDx<-relevel(rse_gene_sACC$PrimaryDx, "Control")


## load
load("rdas/constructed_network_signed_bicor.rda", verbose=TRUE)
net_sACC = net_list[[2]]
net = net_sACC

mergedColors = labels2colors(net$colors)
net$colorsLab = labels2colors(net$colors)
colorDat = data.frame(num = net$colors, col = net$colorsLab, 
	stringsAsFactors=FALSE)
colorDat$Label = paste0("ME", colorDat$num)
colorDat = colorDat[order(colorDat$num),]
colorDat = colorDat[!duplicated(colorDat$num),]
colorDat$numGenes = table(net$colors)[as.character(colorDat$num)]

dim(colorDat)

gList = split(rowData(rse_gene_sACC)$EntrezID, 
	factor(net$colorsLab,levels = colorDat$col))
gList = lapply(gList, function(x) as.character(x[!is.na(x)]))
univ = rowData(rse_gene_sACC)$EntrezID
univ = as.character(univ[!is.na(univ)])

#25 and 4

go_darkgrey = enrichGO(gList$darkgrey, univ = univ,
	OrgDb = "org.Hs.eg.db", ont = "ALL",
	readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 0.05)
go_green = enrichGO(gList$green, univ = univ,
	OrgDb = "org.Hs.eg.db", ont = "ALL",
	readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 0.05)
go_turquoise = enrichGO(gList$turquoise, univ = univ,
	OrgDb = "org.Hs.eg.db", ont = "ALL",
	readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 0.05)
go_magenta = enrichGO(gList$magenta, univ = univ,
	OrgDb = "org.Hs.eg.db", ont = "ALL",
	readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 0.05)


pdf("dotplot/sACC_WGCNA_pathway_enrichment.pdf",h=6,w=12)

dotplot(go_darkgrey, showCategory = 20, title = "sACC darkgrey")
dotplot(go_green, showCategory = 20, title = "sACC green")
dotplot(go_turquoise, showCategory = 20, title = "sACC turquoise")
dotplot(go_magenta, showCategory = 20, title = "sACC magenta")

dev.off()


########################
####Amygdala #####
########################
rm(list = ls())

load(here('exprs_cutoff', 'rse_gene.Rdata'), verbose = TRUE)

rse_gene <- rse_gene[, rse_gene$PrimaryDx %in% c("Control", "MDD")]
rse_gene$Dx <- droplevels(rse_gene$PrimaryDx)
rse_gene_Amygdala <- rse_gene[, rse_gene$BrainRegion %in% c("Amygdala")]
rse_gene_Amygdala$PrimaryDx <- droplevels(rse_gene_Amygdala$PrimaryDx)
rse_gene_Amygdala$PrimaryDx<-relevel(rse_gene_Amygdala$PrimaryDx, "Control")


## load
load("rdas/constructed_network_signed_bicor.rda", verbose=TRUE)
net_Amygdala = net_list[[1]]
net = net_Amygdala

mergedColors = labels2colors(net$colors)
net$colorsLab = labels2colors(net$colors)
colorDat = data.frame(num = net$colors, col = net$colorsLab, 
	stringsAsFactors=FALSE)
colorDat$Label = paste0("ME", colorDat$num)
colorDat = colorDat[order(colorDat$num),]
colorDat = colorDat[!duplicated(colorDat$num),]
colorDat$numGenes = table(net$colors)[as.character(colorDat$num)]

dim(colorDat)
#33 4


gList = split(rowData(rse_gene_Amygdala)$EntrezID, 
	factor(net$colorsLab,levels = colorDat$col))
gList = lapply(gList, function(x) as.character(x[!is.na(x)]))
univ = rowData(rse_gene_Amygdala)$EntrezID
univ = as.character(univ[!is.na(univ)])

go_magenta = enrichGO(gList$magenta, univ = univ,
	OrgDb = "org.Hs.eg.db", ont = "ALL",
	readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 0.05)
go_lightcyan = enrichGO(gList$lightcyan, univ = univ,
	OrgDb = "org.Hs.eg.db", ont = "ALL",
	readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 0.05)
go_darkturquoise = enrichGO(gList$darkturquoise, univ = univ,
	OrgDb = "org.Hs.eg.db", ont = "ALL",
	readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 0.05)

pdf("dotplot/Amygdala_WGCNA_pathway_enrichment.pdf",h=6,w=12)

dotplot(go_magenta, showCategory = 20, title = "Amygdala magenta")
dotplot(go_lightcyan, showCategory = 20, title = "Amygdala lightcyan")
dotplot(go_darkturquoise, showCategory = 20, title = "Amygdala darkturquoise")

dev.off()




#####dotplot code to compare across categories ####### 


c("magenta", "lightcyan", "darkturquoise")
gList = split(rowData(rse_gene_Amygdala)$EntrezID, 
	factor(net$colorsLab,levels = c("magenta", "lightcyan", "darkturquoise")))

go = compareCluster(gList, univ = univ,
	OrgDb = "org.Hs.eg.db", ont = "ALL",
	readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 0.05)

pdf("dotplot/Amygdala_WGCNA_combined_pathway_enrichment.pdf",h=6,w=12)
dotplot(go, showCategory = 20, title = "Amygdala combined")
dev.off()


