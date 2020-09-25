#####################################
# from /dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/R_scripts/DEG_scripts/qSV_model_analysis_DLPFC_PTSD.R
#####################################

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)
library(edgeR)
library(VennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)
library(sessioninfo)

########

load("qSVA_MDD_gene_DEresults.rda", verbose=TRUE)
load("qSVA_MDD_exon_DEresults.rda", verbose=TRUE)
load("qSVA_MDD_jxn_DEresults.rda", verbose=TRUE)
load("qSVA_MDD_tx_DEresults.rda", verbose=TRUE)

# ### save as CSV
# write.csv(outGene_Amyg, file="tables/de_gene_amyg_n541.csv")
# write.csv(outGene_sACC, file="tables/de_gene_sacc_n552.csv")
# write.csv(outExon_Amyg, file="tables/de_exon_amyg_n541.csv")
# write.csv(outExon_sACC, file="tables/de_exon_sacc_n552.csv")
# write.csv(outJxn_Amyg, file="tables/de_jxn_amyg_n541.csv")
# write.csv(outJxn_sACC, file="tables/de_jxn_sacc_n552.csv")
# write.csv(outTx_Amyg, file="tables/de_tx_amyg_n541.csv")
# write.csv(outTx_sACC, file="tables/de_tx_sacc_n552.csv")




## significant rows
sigGene_Amyg_Cnt = outGene_Amyg$gencodeID[outGene_Amyg$q_PrimaryDxControl < 0.05]
sigGene_Amyg_BP = outGene_Amyg$gencodeID[outGene_Amyg$q_PrimaryDxBipolar < 0.05]
sigGene_sACC_Cnt = outGene_sACC$gencodeID[outGene_sACC$q_PrimaryDxControl < 0.05]
sigGene_sACC_BP = outGene_sACC$gencodeID[outGene_sACC$q_PrimaryDxBipolar < 0.05]

sigExon_Amyg_Cnt = outExon_Amyg$gencodeID[outExon_Amyg$q_PrimaryDxControl < 0.05]
sigExon_Amyg_BP = outExon_Amyg$gencodeID[outExon_Amyg$q_PrimaryDxBipolar < 0.05]
sigExon_sACC_Cnt = outExon_sACC$gencodeID[outExon_sACC$q_PrimaryDxControl < 0.05]
sigExon_sACC_BP = outExon_sACC$gencodeID[outExon_sACC$q_PrimaryDxBipolar < 0.05]

# not novel
sigJxn_Amyg_Cnt = outJxn_Amyg_anno$newGeneID[outJxn_Amyg_anno$q_PrimaryDxControl < 0.05]
sigJxn_Amyg_Cnt = sigJxn_Amyg_Cnt[!is.na(sigJxn_Amyg_Cnt)]
sigJxn_Amyg_BP = outJxn_Amyg_anno$newGeneID[outJxn_Amyg_anno$q_PrimaryDxBipolar < 0.05]
sigJxn_Amyg_BP = sigJxn_Amyg_BP[!is.na(sigJxn_Amyg_BP)]
sigJxn_sACC_Cnt = outJxn_sACC_anno$newGeneID[outJxn_sACC_anno$q_PrimaryDxControl < 0.05]
sigJxn_sACC_Cnt = sigJxn_sACC_Cnt[!is.na(sigJxn_sACC_Cnt)]
sigJxn_sACC_BP = outJxn_sACC_anno$newGeneID[outJxn_sACC_anno$q_PrimaryDxBipolar < 0.05]
sigJxn_sACC_BP = sigJxn_sACC_BP[!is.na(sigJxn_sACC_BP)]

sigTx_Amyg_Cnt = outTx_Amyg$gene_id[outTx_Amyg$q_PrimaryDxControl < 0.05]
sigTx_Amyg_BP = outTx_Amyg$gene_id[outTx_Amyg$q_PrimaryDxBipolar < 0.05]
sigTx_sACC_Cnt = outTx_sACC$gene_id[outTx_sACC$q_PrimaryDxControl < 0.05]
sigTx_sACC_BP = outTx_sACC$gene_id[outTx_sACC$q_PrimaryDxBipolar < 0.05]



## gene
venn.diagram(list(Amygdala = sigGene_Amyg_Cnt, sACC = sigGene_sACC_Cnt ), 
             fill = c("slateblue", "skyblue3"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
             margin = .1, imagetype="png",  filename = "venns/venn_fdr05_gene_mddVScnt.png")
venn.diagram(list(Amygdala = sigGene_Amyg_BP, sACC = sigGene_sACC_BP ), 
             fill = c("slateblue", "skyblue3"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
             margin = .1, imagetype="png",  filename = "venns/venn_fdr05_gene_mddVSbip.png")

## exon
venn.diagram(list(Amygdala = sigExon_Amyg_Cnt, sACC = sigExon_sACC_Cnt ), 
             fill = c("slateblue", "skyblue3"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
             margin = .1, imagetype="png",  filename = "venns/venn_fdr05_exon_mddVScnt.png")
venn.diagram(list(Amygdala = sigExon_Amyg_BP, sACC = sigExon_sACC_BP ), 
             fill = c("slateblue", "skyblue3"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
             margin = .1, imagetype="png",  filename = "venns/venn_fdr05_exon_mddVSbip.png")

## jxn
venn.diagram(list(Amygdala = sigJxn_Amyg_Cnt, sACC = sigJxn_sACC_Cnt ), 
             fill = c("slateblue", "skyblue3"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
             margin = .1, imagetype="png",  filename = "venns/venn_fdr05_jxn_mddVScnt.png")
venn.diagram(list(Amygdala = sigJxn_Amyg_BP, sACC = sigJxn_sACC_BP ), 
             fill = c("slateblue", "skyblue3"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
             margin = .1, imagetype="png",  filename = "venns/venn_fdr05_jxn_mddVSbip.png")

## tx
venn.diagram(list(Amygdala = sigTx_Amyg_Cnt, sACC = sigTx_sACC_Cnt ), 
             fill = c("slateblue", "skyblue3"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
             margin = .1, imagetype="png",  filename = "venns/venn_fdr05_tx_mddVScnt.png")
venn.diagram(list(Amygdala = sigTx_Amyg_BP, sACC = sigTx_sACC_BP ), 
             fill = c("slateblue", "skyblue3"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
             margin = .1, imagetype="png",  filename = "venns/venn_fdr05_tx_mddVSbip.png")


## Unique genes of any feature type
Amyg_Cnt = c(sigGene_Amyg_Cnt, sigExon_Amyg_Cnt, sigJxn_Amyg_Cnt, sigTx_Amyg_Cnt)
Amyg_Bp = c(sigGene_Amyg_BP, sigExon_Amyg_BP, sigJxn_Amyg_BP, sigTx_Amyg_BP)

sACC_Cnt = c(sigGene_sACC_Cnt, sigExon_sACC_Cnt, sigJxn_sACC_Cnt, sigTx_sACC_Cnt)
sACC_Bp = c(sigGene_sACC_BP, sigExon_sACC_BP, sigJxn_sACC_BP, sigTx_sACC_BP)

venn.diagram(list(Amygdala = Amyg_Cnt, sACC = sACC_Cnt ),
             fill = c("slateblue", "skyblue3"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
             margin = .1, imagetype="png",  filename = "venns/venn_fdr05_uniquegenes4feat_mddVScnt.png")
venn.diagram(list(Amygdala = Amyg_Bp, sACC = sACC_Bp ), 
             fill = c("slateblue", "skyblue3"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
             margin = .1, imagetype="png",  filename = "venns/venn_fdr05_uniquegenes4feat_mddVSbip.png")



#### run enrichment analysis ####

moduleGeneList = outGene_sACC$EntrezID[outGene_sACC$q_PrimaryDxControl<0.05]
moduleGeneList = moduleGeneList[!is.na(moduleGeneList)]
length(moduleGeneList)
# 588

geneUniverse = as.character(outGene_sACC$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

## run enrichment analysis
goBP <- enrichGO(moduleGeneList,
                 universe = geneUniverse, OrgDb = org.Hs.eg.db,
                 ont = "BP", pAdjustMethod = "BH",
                 pvalueCutoff  = .21, qvalueCutoff  = .75,
                 readable= TRUE)
goMF <- enrichGO(moduleGeneList,
                 universe = geneUniverse, OrgDb = org.Hs.eg.db,
                 ont = "MF", pAdjustMethod = "BH",
                 pvalueCutoff  = .21, qvalueCutoff  = .75,
                 readable= TRUE)
goCC <- enrichGO(moduleGeneList,
                 universe = geneUniverse, OrgDb = org.Hs.eg.db,
                 ont = "CC", pAdjustMethod = "BH",
                 pvalueCutoff  = .21, qvalueCutoff  = .75,
                 readable= TRUE)

pdf("gene_enrichments_sACC_fdr05.pdf",h=6,w=12)
dotplot(goBP, title="Biological Processes, sACC")
dotplot(goMF, title="Molecular Functions, sACC")
dotplot(goCC, title="Cellular Components, sACC")
dev.off()


### split direction (use FDR 0.1)
up = outGene_sACC$EntrezID[outGene_sACC$q_PrimaryDxControl<0.1 & outGene_sACC$PrimaryDxBipolar>0]
down = outGene_sACC$EntrezID[outGene_sACC$q_PrimaryDxControl<0.1 & outGene_sACC$PrimaryDxBipolar<0]
moduleGeneList = list(up, down)
moduleGeneList = lapply(moduleGeneList, function(x) x[!is.na(x)])
names(moduleGeneList) = c("CNT>MDD","CNT<MDD")

lapply(moduleGeneList, length)
# 509, 555

geneUniverse = as.character(outGene_sACC$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

## run enrichment analysis
goBP <- enrichGO(moduleGeneList[[1]],
                 universe = geneUniverse, OrgDb = org.Hs.eg.db,
                 ont = "BP", pAdjustMethod = "BH",
                 pvalueCutoff  = .2, qvalueCutoff  = .5,
                 readable= TRUE)
goBP2 <- enrichGO(moduleGeneList[[2]],
                  universe = geneUniverse, OrgDb = org.Hs.eg.db,
                  ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff  = .2, qvalueCutoff  = .5,
                  readable= TRUE)
goMF <- enrichGO(moduleGeneList[[1]],
                 universe = geneUniverse, OrgDb = org.Hs.eg.db,
                 ont = "MF", pAdjustMethod = "BH",
                 pvalueCutoff  = .2, qvalueCutoff  = .5,
                 readable= TRUE)
goMF2 <- enrichGO(moduleGeneList[[2]],
                  universe = geneUniverse, OrgDb = org.Hs.eg.db,
                  ont = "MF", pAdjustMethod = "BH",
                  pvalueCutoff  = .2, qvalueCutoff  = .5,
                  readable= TRUE)
goCC <- enrichGO(moduleGeneList[[1]],
                 universe = geneUniverse, OrgDb = org.Hs.eg.db,
                 ont = "CC", pAdjustMethod = "BH",
                 pvalueCutoff  = .2, qvalueCutoff  = .5,
                 readable= TRUE)
goCC2 <- enrichGO(moduleGeneList[[2]],
                  universe = geneUniverse, OrgDb = org.Hs.eg.db,
                  ont = "CC", pAdjustMethod = "BH",
                  pvalueCutoff  = .2, qvalueCutoff  = .5,
                  readable= TRUE)

pdf("gene_enrichments_sACC_up_down_fdr10.pdf",h=6,w=12)
dotplot(goBP, title="CNT>MDD, Biological Processes, sACC")
dotplot(goBP2, title="CNT<MDD, Biological Processes, sACC")
dotplot(goMF, title="CNT>MDD, Molecular Functions, sACC")
dotplot(goMF2, title="CNT<MDD, Molecular Functions, sACC")
dotplot(goCC, title="CNT>MDD, Cellular Components, sACC")
dotplot(goCC2, title="CNT<MDD, Cellular Components, sACC")
dev.off()




#### Amygdala

moduleGeneList = outGene_Amyg$EntrezID[outGene_Amyg$q_PrimaryDxControl<0.05]
moduleGeneList = moduleGeneList[!is.na(moduleGeneList)]
length(moduleGeneList)
# 77

geneUniverse = as.character(outGene_Amyg$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

## run enrichment analysis
goBP <- enrichGO(moduleGeneList,
                 universe = geneUniverse, OrgDb = org.Hs.eg.db,
                 ont = "BP", pAdjustMethod = "BH",
                 pvalueCutoff  = .21, qvalueCutoff  = .75,
                 readable= TRUE)
goMF <- enrichGO(moduleGeneList,
                 universe = geneUniverse, OrgDb = org.Hs.eg.db,
                 ont = "MF", pAdjustMethod = "BH",
                 pvalueCutoff  = .21, qvalueCutoff  = .75,
                 readable= TRUE)
goCC <- enrichGO(moduleGeneList,
                 universe = geneUniverse, OrgDb = org.Hs.eg.db,
                 ont = "CC", pAdjustMethod = "BH",
                 pvalueCutoff  = .21, qvalueCutoff  = .75,
                 readable= TRUE)

pdf("gene_enrichments_Amyg_fdr05.pdf",h=6,w=12)
dotplot(goBP, title="Biological Processes, Amygdala")
dotplot(goMF, title="Molecular Functions, Amygdala")
dotplot(goCC, title="Cellular Components, Amygdala")
dev.off()


## split up down
up = outGene_Amyg$EntrezID[outGene_Amyg$q_PrimaryDxControl<0.1 & outGene_Amyg$PrimaryDxBipolar>0]
down = outGene_Amyg$EntrezID[outGene_Amyg$q_PrimaryDxControl<0.1 & outGene_Amyg$PrimaryDxBipolar<0]
moduleGeneList = list(up, down)
moduleGeneList = lapply(moduleGeneList, function(x) x[!is.na(x)])
names(moduleGeneList) = c("CNT>MDD","CNT<MDD")

lapply(moduleGeneList, length)
# 108, 64

geneUniverse = as.character(outGene_Amyg$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

## run enrichment analysis
goBP <- enrichGO(moduleGeneList[[1]],
                 universe = geneUniverse, OrgDb = org.Hs.eg.db,
                 ont = "BP", pAdjustMethod = "BH",
                 pvalueCutoff  = .2, qvalueCutoff  = .5,
                 readable= TRUE)
goBP2 <- enrichGO(moduleGeneList[[2]],
                  universe = geneUniverse, OrgDb = org.Hs.eg.db,
                  ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff  = .2, qvalueCutoff  = .5,
                  readable= TRUE)
goMF <- enrichGO(moduleGeneList[[1]],
                 universe = geneUniverse, OrgDb = org.Hs.eg.db,
                 ont = "MF", pAdjustMethod = "BH",
                 pvalueCutoff  = .2, qvalueCutoff  = .5,
                 readable= TRUE)
goMF2 <- enrichGO(moduleGeneList[[2]],
                  universe = geneUniverse, OrgDb = org.Hs.eg.db,
                  ont = "MF", pAdjustMethod = "BH",
                  pvalueCutoff  = .2, qvalueCutoff  = .5,
                  readable= TRUE)
goCC <- enrichGO(moduleGeneList[[1]],
                 universe = geneUniverse, OrgDb = org.Hs.eg.db,
                 ont = "CC", pAdjustMethod = "BH",
                 pvalueCutoff  = .2, qvalueCutoff  = .5,
                 readable= TRUE)
goCC2 <- enrichGO(moduleGeneList[[2]],
                  universe = geneUniverse, OrgDb = org.Hs.eg.db,
                  ont = "CC", pAdjustMethod = "BH",
                  pvalueCutoff  = .2, qvalueCutoff  = .5,
                  readable= TRUE)

pdf("gene_enrichments_Amyg_up_down_fdr10.pdf",h=6,w=12)
dotplot(goBP, title="CNT>MDD, Biological Processes, Amygdala")
dotplot(goBP2, title="CNT<MDD, Biological Processes, Amygdala")
dotplot(goMF, title="CNT>MDD, Molecular Functions, Amygdala")
dotplot(goMF2, title="CNT<MDD, Molecular Functions, Amygdala")
dotplot(goCC, title="CNT>MDD, Cellular Components, Amygdala")
dotplot(goCC2, title="CNT<MDD, Cellular Components, Amygdala")
dev.off()


#sgejobs::job_single('qSV_model_DE_explore', create_shell = TRUE, memory = '80G', command = "Rscript qSV_model_DE_explore.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
