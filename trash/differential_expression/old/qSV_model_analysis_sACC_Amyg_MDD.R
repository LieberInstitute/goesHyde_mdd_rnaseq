#####################################
# from /dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/R_scripts/DEG_scripts/qSV_model_analysis_DLPFC_PTSD.R
#####################################

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)
library(edgeR)

setwd('/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/differential_expression/')

###############################
##### Load, clean, combine ####
###############################

#load objects
load('../preprocessed_data/rse_gene_goesHyde_MDD_n634.Rdata', verbose=TRUE)
rse_mdd = rse_gene
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseGene_n511.rda", verbose=TRUE)
rse_bip = rse_gene

## combine
# make colData consistent
rse_mdd$AgeDeath = rse_mdd$Age
rse_mdd$RNum = rse_mdd$SAMPLE_ID
rse_mdd$BrNum = as.character(rse_mdd$BrNum)
colKeep = c("RNum","BrNum","Sex","Race","AgeDeath","BrainRegion","RIN","PrimaryDx","overallMapRate","totalAssignedGene","mitoRate","rRNA_rate")
colData(rse_mdd) = colData(rse_mdd)[,colKeep]
colData(rse_bip) = colData(rse_bip)[,colKeep]

# make rowData consistent
rowData(rse_bip)$Symbol = rowData(rse_mdd)$Symbol 	# fill in blank symbols
rowData(rse_bip)$meanExprs = rowData(rse_bip)$gencodeTx = NULL
rowData(rse_mdd)$meanExprs = rowData(rse_mdd)$gencodeTx = NULL

### combine
rse_both = cbind(rse_mdd, rse_bip)

## drop
qc = read.csv("../qc_checks/qc_dropping_results.csv", stringsAsFactors = FALSE)
qc = qc[rowSums(qc[,13:16])>0,]
rse_both = rse_both[,-which(rse_both$RNum %in% qc$SAMPLE_ID | rse_both$RNum %in% c("R17538","R18853") | 
						rse_both$PrimaryDx == "Other" | 
						rse_both$overallMapRate <0.5) ]
rse_both$PrimaryDx = droplevels(rse_both$PrimaryDx)
rse_both$PrimaryDx = relevel(rse_both$PrimaryDx, ref="MDD")
pd = colData(rse_both)
table(pd$BrainRegion,pd$PrimaryDx)
           # MDD Control Bipolar
  # Amygdala 236     190     120
  # sACC     228     204     125

table(pd$Sex, pd$PrimaryDx)
    # MDD Control Bipolar
  # F 155      78      96
  # M 309     316     149

summary(pd$AgeDeath)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 17.37   34.27   47.21   46.50   55.86   95.27 


#load cov_rse object
load("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/data/degradation_rse_MDDseq_BothRegions.rda")
cov_rsemdd = cov_rse
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/data/degradation_rse_BipSeq_BothRegions.rda")
cov_rsebip = cov_rse

# combine
cov_rsemdd$AgeDeath = cov_rsemdd$Age
cov_rsemdd$RNum = cov_rsemdd$SAMPLE_ID
cov_rsemdd$BrNum = as.character(cov_rsemdd$BrNum)
colKeep = c("RNum","BrNum","Sex","Race","AgeDeath","BrainRegion","RIN","PrimaryDx","overallMapRate","totalAssignedGene","mitoRate","rRNA_rate")
colData(cov_rsemdd) = colData(cov_rsemdd)[,colKeep]
colData(cov_rsebip) = colData(cov_rsebip)[,colKeep]

cov_rseboth = cbind(cov_rsemdd, cov_rsebip)
cov_rseboth = cov_rseboth[,cov_rseboth$RNum %in% rse_both$RNum]

stopifnot(identical(cov_rseboth$RNum, rse_both$RNum))

# objects:
rse_gene = rse_both
cov_rse = cov_rseboth


###########
# filter ##
## gene
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')
geneIndex = rowMeans(assays(rse_gene)$rpkm) > 0.25  ## both regions
rse_gene = rse_gene[geneIndex,]


###################
##### get qSVs
###################

modJoint = model.matrix(~PrimaryDx + AgeDeath + Sex + mitoRate + rRNA_rate + 
				totalAssignedGene + RIN, data = colData(rse_gene))

degExprs = log2(assays(cov_rse)$count+1)
k = num.sv(degExprs, modJoint) # 22
qSV_mat = prcomp(t(degExprs))$x[,1:k]
varExplQsva = getPcaVars(prcomp(t(degExprs)))
varExplQsva[1:k]
sum(varExplQsva[1:k]) # 88.7%

# model w/o interaction to subset by region
modSep = model.matrix(~PrimaryDx + AgeDeath + Sex + mitoRate + rRNA_rate +
				totalAssignedGene + RIN, data=colData(rse_gene))

#########################
## split back by region #
#########################

#### both ####
sACC_Index = which(colData(rse_gene)$BrainRegion == "sACC")
mod_sACC = cbind(modSep[sACC_Index,], qSV_mat[sACC_Index, ])
Amyg_Index = which(colData(rse_gene)$BrainRegion == "Amygdala")
mod_Amyg = cbind(modSep[Amyg_Index,], qSV_mat[Amyg_Index, ])



#################
### Gene ########
#################

##### sACC ######
dge_sACC = DGEList(counts = assays(rse_gene[,sACC_Index])$counts, 
	genes = rowData(rse_gene))
dge_sACC = calcNormFactors(dge_sACC)
vGene_sACC = voom(dge_sACC,mod_sACC, plot=FALSE)

fitGene_sACC = lmFit(vGene_sACC)
eBGene_sACC = eBayes(fitGene_sACC)
outGene_sACC = topTable(eBGene_sACC,coef=2:3,
	p.value = 1,number=nrow(rse_gene))
outGene_sACC = outGene_sACC[rownames(rse_gene),]
	
## significance levels
pvalMat = as.matrix(eBGene_sACC$p.value)[,2:3]
qvalMat = pvalMat
qvalMat[,1:2] = p.adjust(pvalMat[,1:2],method="fdr") 
colnames(pvalMat) = paste0("P_",colnames(pvalMat))
colnames(qvalMat) = paste0("q_",colnames(qvalMat))

outGene_sACC = cbind(outGene_sACC,cbind(pvalMat, qvalMat))
sum(outGene_sACC$q_PrimaryDxControl < 0.05)
# 739

# sum(outGene_sACC$q_PrimaryDxControl < 0.05)
# [1] 739
# sum(outGene_sACC$q_PrimaryDxControl < 0.01)
# [1] 255

# sum(outGene_sACC$q_PrimaryDxBipolar < 0.05)
# [1] 305
# sum(outGene_sACC$q_PrimaryDxBipolar < 0.01)


##### Amygdala ######
dge_Amyg = DGEList(counts = assays(rse_gene[,Amyg_Index])$counts, 
	genes = rowData(rse_gene))
dge_Amyg = calcNormFactors(dge_Amyg)
vGene_Amyg = voom(dge_Amyg,mod_Amyg, plot=FALSE)

fitGene_Amyg = lmFit(vGene_Amyg)
eBGene_Amyg = eBayes(fitGene_Amyg)
outGene_Amyg = topTable(eBGene_Amyg,coef=2:3,
	p.value = 1,number=nrow(rse_gene))
outGene_Amyg = outGene_Amyg[rownames(rse_gene),]

## significance levels
pvalMat = as.matrix(eBGene_Amyg$p.value)[,2:3]
qvalMat = pvalMat
qvalMat[,1:2] = p.adjust(pvalMat[,1:2],method="fdr") 
colnames(pvalMat) = paste0("P_",colnames(pvalMat))
colnames(qvalMat) = paste0("q_",colnames(qvalMat))

outGene_Amyg = cbind(outGene_Amyg,cbind(pvalMat, qvalMat))
sum(outGene_Amyg$q_PrimaryDxControl < 0.05)
# 128

# sum(outGene_Amyg$q_PrimaryDxControl < 0.05)
# [1] 128
# sum(outGene_Amyg$q_PrimaryDxControl < 0.01)
# [1] 30

# sum(outGene_Amyg$q_PrimaryDxBipolar < 0.05)
# [1] 120
# sum(outGene_Amyg$q_PrimaryDxBipolar < 0.01)
# [1] 39



library(VennDiagram)

venn.diagram(list(Amygdala = outGene_Amyg$ensemblID[outGene_Amyg$q_PrimaryDxControl < 0.05], 
				sACC = outGene_sACC$ensemblID[outGene_sACC$q_PrimaryDxControl < 0.05] ), 
	fill = c("slateblue", "skyblue3"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_fdr05_mddVScnt.png")

venn.diagram(list(Amygdala = outGene_Amyg$ensemblID[outGene_Amyg$q_PrimaryDxBipolar < 0.05], 
				sACC = outGene_sACC$ensemblID[outGene_sACC$q_PrimaryDxBipolar < 0.05] ), 
	fill = c("slateblue", "skyblue3"), main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_fdr05_mddVSbip.png")


	
	
orderedSacc = outGene_sACC[order(outGene_sACC$P_PrimaryDxControl),]
orderedAmyg = outGene_Amyg[order(outGene_Amyg$P_PrimaryDxControl),]
write.csv(orderedSacc[which(orderedSacc$q_PrimaryDxControl < 0.05),], file="de_results_sACC.csv")
write.csv(orderedAmyg[which(orderedAmyg$q_PrimaryDxControl < 0.05),], file="de_results_Amyg.csv")


######################
### plot
library(RColorBrewer)

geneRpkm = assays(rse_gene)$rpkm
yExprs = as.matrix(log2(geneRpkm+1))

###### sACC
## residualize expression			
saccExprs = cleaningY(yExprs[,sACC_Index], mod_sACC, P=3)
pdSacc = pd[sACC_Index,]

sigOrderMat = as.data.frame(apply(outGene_sACC[,c(15:16)], 2, 
                                  function(x) order(x)[1:100]))
ooL = sigOrderMat$"P_PrimaryDxControl"

pdf("top_genes_MDD_vs_cnt_sACC.pdf",h=6,w=6)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(4,"Dark2"))
for(i in ooL) {
  boxplot(saccExprs[i,] ~ as.numeric(pdSacc$PrimaryDx) , outline=FALSE, xaxt="n",
	   pch = 21,
       cex=2, xlab="",
       ylab="Residualized Expression",
       main = paste0(outGene_sACC$Symbol[i], "\n", outGene_sACC$gencodeID[i]) )
 points(saccExprs[i,] ~ jitter(as.numeric(factor(pdSacc$PrimaryDx))),
       pch=21, bg=as.numeric(pdSacc$PrimaryDx), cex=1.5)
 axis(1, at=1:1:length(levels(pdSacc$PrimaryDx)), labels = levels(as.factor(pdSacc$PrimaryDx)))
 # abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
 # legend("top", paste0("p=",signif(voomStatsC$"pvalue_NPC-ACC_DORSAL"[i],3)), bg="white")
 # if (i==ooL[1]) { legend("topright", levels(factor(pd$LINE)),
	# pch = 15, col = pal, cex=1) }
}
dev.off()



###### Amygdala
## residualize expression			
amygExprs = cleaningY(yExprs[,Amyg_Index], mod_Amyg, P=3)
pdAmyg = pd[Amyg_Index,]


sigOrderMat = as.data.frame(apply(outGene_Amyg[,c(15:16)], 2, 
                                  function(x) order(x)[1:100]))
ooL = sigOrderMat$"P_PrimaryDxControl"

pdf("top_genes_MDD_vs_cnt_Amygdala.pdf",h=6,w=6)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(4,"Dark2"))
for(i in ooL) {
  boxplot(amygExprs[i,] ~ as.numeric(pdAmyg$PrimaryDx) , outline=FALSE, xaxt="n",
	   pch = 21,
       cex=2, xlab="",
       ylab="Residualized Expression",
       main = paste0(outGene_Amyg$Symbol[i], "\n", outGene_Amyg$gencodeID[i]) )
 points(amygExprs[i,] ~ jitter(as.numeric(factor(pdAmyg$PrimaryDx))),
       pch=21, bg=as.numeric(pdAmyg$PrimaryDx), cex=1.5)
 axis(1, at=1:1:length(levels(pdAmyg$PrimaryDx)), labels = levels(as.factor(pdAmyg$PrimaryDx)))
 # abline(v=c(4.5,5.5,6.5), col="grey", lty=2)
 # legend("top", paste0("p=",signif(voomStatsC$"pvalue_NPC-ACC_DORSAL"[i],3)), bg="white")
 # if (i==ooL[1]) { legend("topright", levels(factor(pd$LINE)),
	# pch = 15, col = pal, cex=1) }
}
dev.off()





	

############################## 
## run enrichment analysis ###
##############################
library(clusterProfiler)
library(org.Hs.eg.db)

sigQsva = outGene_sACC[outGene_sACC$q_PrimaryDxControl< 0.05,]
sigQsvaGenes = unique(as.character(sigQsva$EntrezID[!is.na(sigQsva$EntrezID)]))
length(sigQsvaGenes)
# 665

moduleGeneList = sigQsvaGenes

geneUniverse = as.character(outGene_sACC$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

## run enrichment analysis
goBP <- enrichGO(moduleGeneList,
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "BP", pAdjustMethod = "BH",
                pvalueCutoff  = .1, qvalueCutoff  = .5,
				readable= TRUE)
goMF <- enrichGO(moduleGeneList,
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "MF", pAdjustMethod = "BH",
                pvalueCutoff  = .1, qvalueCutoff  = .5,
				readable= TRUE)
goCC <- enrichGO(moduleGeneList,
                universe = geneUniverse, OrgDb = org.Hs.eg.db,
                ont = "CC", pAdjustMethod = "BH",
                pvalueCutoff  = .1, qvalueCutoff  = .5,
				readable= TRUE)
kegg <- enrichGO(moduleGeneList,
                universe = geneUniverse,  pAdjustMethod = "BH",
                pvalueCutoff  = .1, qvalueCutoff  = .5)

pdf("de_up_down_enrichments.pdf",h=5,w=9)
dotplot(goBP, includeAll="TRUE", title="Biological Processes")
dotplot(goMF, includeAll="TRUE", title="Molecular Functions")
dotplot(goCC, includeAll="TRUE", title="Cellular Components")
dotplot(kegg, includeAll="TRUE", title="KEGG terms")
dev.off()









## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
