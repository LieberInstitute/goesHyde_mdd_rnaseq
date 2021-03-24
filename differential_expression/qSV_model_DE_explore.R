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
library(tidyverse)
library(GGally)
library(viridis)
library(plotly)
library(here)

#### Load Data ####
data_type <- c("gene","exon","jxn","tx")
names(data_type) <- data_type

paths <- map(data_type, ~here("differential_expression","data",paste0("qSVA_MDD_", .x, "_DEresults.rda")))
load(paths$gene, verbose = TRUE)
# outGene <- list(amyg = outGene_Amyg, sacc = outGene_sACC)

#### Find Significant rows ####
dx_test <- c(ctrl = "q_PrimaryDxControl", bp = "q_PrimaryDxBipolar")

get_sig <- function(outdf, cols = dx_test , cutoff = 0.05){
  sig <- map(cols, ~outdf$gencodeID[outdf[[.x]] < cutoff])
  return(sig)
}

sigGene <- map(outGene, get_sig)
map(sigGene, ~map_int(.x,length))

sigGene_ilr <- map(outGene_ilr, get_sig)
map(sigGene_ilr, ~map(.x,length))


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


#### venn diagrams ####

de_venn <- function(in_list, fn){
  full_fn = here("differential_expression","venns",paste0(fn,".png"))


  venn.diagram(in_list,
               fill = c("slateblue", "skyblue3"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
               margin = .1, imagetype="png",  filename = full_fn)
}

de_venn(map(sigGene, ~pluck(.x, "ctrl")), fn = "venn_fdr05_gene_mddVSctrl")
de_venn(map(sigGene, ~pluck(.x, "bp")), fn = "venn_fdr05_gene_mddVSbip")

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
# 514

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

map_int(c(goBP,goMF, goCC), nrow)
# [1] 0 0 6

pdf("gene_enrichments_sACC_fdr05.pdf",h=6,w=12)
dotplot(goBP, title="Biological Processes, sACC")
dotplot(goMF, title="Molecular Functions, sACC")
dotplot(goCC, title="Cellular Components, sACC")
dev.off()
save(goBP, goMF, goCC, file = "gene_enrichments_sACC_fdr05.Rdata")


### split direction (use FDR 0.1)
up = outGene_sACC$EntrezID[outGene_sACC$q_PrimaryDxControl<0.1 & outGene_sACC$PrimaryDxBipolar>0]
down = outGene_sACC$EntrezID[outGene_sACC$q_PrimaryDxControl<0.1 & outGene_sACC$PrimaryDxBipolar<0]
moduleGeneList = list(up, down)
moduleGeneList = lapply(moduleGeneList, function(x) x[!is.na(x)])
names(moduleGeneList) = c("CNT>MDD","CNT<MDD")

map_int(moduleGeneList, length)
# CNT>MDD CNT<MDD
# 494     497

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

map_int(c(goBP, goBP2, goMF, goMF2, goCC, goCC2), nrow)
# [1] 0 0 3 0 0 0

pdf("gene_enrichments_sACC_up_down_fdr10.pdf",h=6,w=12)
dotplot(goBP, title="CNT>MDD, Biological Processes, sACC")
dotplot(goBP2, title="CNT<MDD, Biological Processes, sACC")
dotplot(goMF, title="CNT>MDD, Molecular Functions, sACC")
dotplot(goMF2, title="CNT<MDD, Molecular Functions, sACC")
dotplot(goCC, title="CNT>MDD, Cellular Components, sACC")
dotplot(goCC2, title="CNT<MDD, Cellular Components, sACC")
dev.off()
save(goBP, goBP2, goMF, goMF2, goCC, goCC2,file="gene_enrichments_sACC_up_down_fdr10.Rdata")



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

map_int(c(goBP,goMF, goCC), nrow)
# [1] 93 27 28

pdf("gene_enrichments_Amyg_fdr05.pdf",h=6,w=12)
dotplot(goBP, title="Biological Processes, Amygdala")
dotplot(goMF, title="Molecular Functions, Amygdala")
dotplot(goCC, title="Cellular Components, Amygdala")
dev.off()
save(goBP, goMF, goCC, file = "gene_enrichments_Amyg_fdr05.Rdata")

## split up down
up = outGene_Amyg$EntrezID[outGene_Amyg$q_PrimaryDxControl<0.1 & outGene_Amyg$PrimaryDxBipolar>0]
down = outGene_Amyg$EntrezID[outGene_Amyg$q_PrimaryDxControl<0.1 & outGene_Amyg$PrimaryDxBipolar<0]
moduleGeneList = list(up, down)
moduleGeneList = lapply(moduleGeneList, function(x) x[!is.na(x)])
names(moduleGeneList) = c("CNT>MDD","CNT<MDD")

map_int(moduleGeneList, length)
# CNT>MDD CNT<MDD
# 121      71

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

map_int(c(goBP, goBP2, goMF, goMF2, goCC, goCC2), nrow)
# [1] 200  21  45   0  28   0

pdf("gene_enrichments_Amyg_up_down_fdr10.pdf",h=6,w=12)
dotplot(goBP, title="CNT>MDD, Biological Processes, Amygdala")
dotplot(goBP2, title="CNT<MDD, Biological Processes, Amygdala")
dotplot(goMF, title="CNT>MDD, Molecular Functions, Amygdala")
dotplot(goMF2, title="CNT<MDD, Molecular Functions, Amygdala")
dotplot(goCC, title="CNT>MDD, Cellular Components, Amygdala")
dotplot(goCC2, title="CNT<MDD, Cellular Components, Amygdala")
dev.off()
save(goBP, goBP2, goMF, goMF2, goCC, goCC2, file = "gene_enrichments_Amyg_up_down_fdr10.Rdata")


#sgejobs::job_single('qSV_model_DE_explore', create_shell = TRUE, memory = '80G', command = "Rscript qSV_model_DE_explore.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
