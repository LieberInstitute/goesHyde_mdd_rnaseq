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

#load objects
load('../exprs_cutoff/rse_gene.Rdata', verbose=TRUE)
load('../exprs_cutoff/rse_exon.Rdata', verbose=TRUE)
load('../exprs_cutoff/rse_jxn.Rdata', verbose=TRUE)
load('../exprs_cutoff/rse_tx.Rdata', verbose=TRUE)


##############################
############ load degredation

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
cov_rseboth = cov_rseboth[,match(rse_gene$RNum,cov_rseboth$RNum)]	# put in order / drop dropped samples
stopifnot(identical(cov_rseboth$RNum, rse_gene$RNum))

cov_rse = cov_rseboth


###################
##### get qSVs
###################

modJoint = model.matrix(~PrimaryDx + AgeDeath + Sex + mitoRate + rRNA_rate + 
				totalAssignedGene + RIN, data = colData(rse_gene))

degExprs = log2(assays(cov_rse)$count+1)
k = num.sv(degExprs, modJoint) # 23
qSV_mat = prcomp(t(degExprs))$x[,1:k]
varExplQsva = getPcaVars(prcomp(t(degExprs)))
varExplQsva[1:k]
sum(varExplQsva[1:k]) # 87.888%

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



###################
## load DE results
########

load("qSVA_MDD_gene_DEresults.rda", verbose=TRUE)
load("qSVA_MDD_exon_DEresults.rda", verbose=TRUE)
load("qSVA_MDD_jxn_DEresults.rda", verbose=TRUE)
load("qSVA_MDD_tx_DEresults.rda", verbose=TRUE)




##################
## expression

rse_gene$Dx = as.factor(rse_gene$PrimaryDx)
rse_gene$Dx = factor(rse_gene$Dx, levels(rse_gene$Dx)[c(2,3,1)] )

pdSacc = colData(rse_gene)[sACC_Index,]
pdAmyg = colData(rse_gene)[Amyg_Index,]

gRpkm = recount::getRPKM(rse_gene,"Length")
eRpkm = recount::getRPKM(rse_exon,"Length")
jRpkm = recount::getRPKM(rse_jxn,"Length")
tRpkm = assays(rse_tx)$tpm

gExprs = as.matrix(log2(gRpkm+1))
eExprs = as.matrix(log2(eRpkm+1))
jExprs = as.matrix(log2(jRpkm+1))
tExprs = as.matrix(log2(tRpkm+1))

gSaccExprs = cleaningY(gExprs[,sACC_Index], mod_sACC, P=3)
gAmygExprs = cleaningY(gExprs[,Amyg_Index], mod_Amyg, P=3)

eSaccExprs = cleaningY(eExprs[,sACC_Index], mod_sACC, P=3)
eAmygExprs = cleaningY(eExprs[,Amyg_Index], mod_Amyg, P=3)

jSaccExprs = cleaningY(jExprs[,sACC_Index], mod_sACC, P=3)
jAmygExprs = cleaningY(jExprs[,Amyg_Index], mod_Amyg, P=3)

tSaccExprs = cleaningY(tExprs[,sACC_Index], mod_sACC, P=3)
tAmygExprs = cleaningY(tExprs[,Amyg_Index], mod_Amyg, P=3)






######################
### plot
library(RColorBrewer)


###### sACC
## residualize expression			

sigOrderMat = as.data.frame(apply(outGene_sACC[,c(17:18)], 2, function(x) order(x)[1:100]))
ooL = sigOrderMat$"P_PrimaryDxControl"
pdf("plots/top_genes_MDD_vs_cnt_sACC.pdf",h=6,w=6)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(4,"Dark2"))
for(i in ooL) {
  boxplot(gSaccExprs[i,] ~ as.numeric(pdSacc$Dx) , outline=FALSE, xaxt="n",
	   pch = 21,
       cex=2, xlab="",
	   ylim = c(min(gSaccExprs[i,]),max(gSaccExprs[i,])),
       ylab="Residualized Expression",
       main = paste0(outGene_sACC$Symbol[i], "\n", outGene_sACC$gencodeID[i]) )
 points(gSaccExprs[i,] ~ jitter(as.numeric(factor(pdSacc$Dx))),
       pch=21, bg=as.numeric(pdSacc$Dx), cex=1.5)
 axis(1, at=1:1:length(levels(pdSacc$Dx)), labels = levels(as.factor(pdSacc$Dx)))
 legend("top", paste0("p=",signif(outGene_sACC$P_PrimaryDxControl[i],3)), bg="white")
}
dev.off()


sigOrderMat = as.data.frame(apply(outExon_sACC[,c(23:24)], 2, function(x) order(x)[1:100]))
ooL = sigOrderMat$"P_PrimaryDxControl"
pdf("plots/top_exons_MDD_vs_cnt_sACC.pdf",h=6,w=6)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(4,"Dark2"))
for(i in ooL) {
  boxplot(eSaccExprs[i,] ~ as.numeric(pdSacc$Dx) , outline=FALSE, xaxt="n",
	   pch = 21,
       cex=2, xlab="",
	   ylim = c(min(eSaccExprs[i,]),max(eSaccExprs[i,])),
       ylab="Residualized Expression",
       main = paste0(outExon_sACC$Symbol[i], "\n", outExon_sACC$coord[i]) )
 points(eSaccExprs[i,] ~ jitter(as.numeric(factor(pdSacc$Dx))),
       pch=21, bg=as.numeric(pdSacc$Dx), cex=1.5)
 axis(1, at=1:1:length(levels(pdSacc$Dx)), labels = levels(as.factor(pdSacc$Dx)))
 legend("top", paste0("p=",signif(outExon_sACC$P_PrimaryDxControl[i],3)), bg="white")
}
dev.off()


sigOrderMat = as.data.frame(apply(outJxn_sACC_anno[,c(25:26)], 2, function(x) order(x)[1:100]))
ooL = sigOrderMat$"P_PrimaryDxControl"
pdf("plots/top_jxns_MDD_vs_cnt_sACC.pdf",h=6,w=6)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(4,"Dark2"))
for(i in ooL) {
  jxnI = rownames(outJxn_sACC_anno)[i]
  indI = which(rownames(jSaccExprs) == jxnI)
  boxplot(jSaccExprs[indI,] ~ as.numeric(pdSacc$Dx) , outline=FALSE, xaxt="n",
	   pch = 21,
       cex=2, xlab="",
	   ylim = c(min(jSaccExprs[indI,]),max(jSaccExprs[indI,])),
       ylab="Residualized Expression",
       main = paste0(outJxn_sACC_anno$newGeneSymbol[i],", ",outJxn_sACC_anno$Class[i],"\n", rownames(outJxn_sACC_anno)[i]) )
 points(jSaccExprs[indI,] ~ jitter(as.numeric(factor(pdSacc$Dx))),
       pch=21, bg=as.numeric(pdSacc$Dx), cex=1.5)
 axis(1, at=1:1:length(levels(pdSacc$Dx)), labels = levels(as.factor(pdSacc$Dx)))
 legend("top", paste0("p=",signif(outJxn_sACC_anno$P_PrimaryDxControl[i],3)), bg="white")
}
dev.off()


sigOrderMat = as.data.frame(apply(outTx_sACC[,c(11:12)], 2, function(x) order(x)[1:100]))
ooL = sigOrderMat$"P_PrimaryDxControl"
pdf("plots/top_tx_MDD_vs_cnt_sACC.pdf",h=6,w=6)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(4,"Dark2"))
for(i in ooL) {
  boxplot(tSaccExprs[i,] ~ as.numeric(pdSacc$Dx) , outline=FALSE, xaxt="n",
	   pch = 21,
       cex=2, xlab="",
	   ylim = c(min(tSaccExprs[i,]),max(tSaccExprs[i,])),
       ylab="Residualized Expression",
       main = paste0(outTx_sACC$gene_name[i], "\n", rownames(outTx_sACC)[i]) )
 points(tSaccExprs[i,] ~ jitter(as.numeric(factor(pdSacc$Dx))),
       pch=21, bg=as.numeric(pdSacc$Dx), cex=1.5)
 axis(1, at=1:1:length(levels(pdSacc$Dx)), labels = levels(as.factor(pdSacc$Dx)))
 legend("top", paste0("p=",signif(outTx_sACC$P_PrimaryDxControl[i],3)), bg="white")
}
dev.off()








###### Amygdala
## residualize expression			

sigOrderMat = as.data.frame(apply(outGene_Amyg[,c(17:18)], 2, function(x) order(x)[1:100]))
ooL = sigOrderMat$"P_PrimaryDxControl"
pdf("plots/top_genes_MDD_vs_cnt_Amyg.pdf",h=6,w=6)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(4,"Dark2"))
for(i in ooL) {
  boxplot(gAmygExprs[i,] ~ as.numeric(pdAmyg$Dx) , outline=FALSE, xaxt="n",
	   pch = 21,
       cex=2, xlab="",
	   ylim = c(min(gAmygExprs[i,]),max(gAmygExprs[i,])),
       ylab="Residualized Expression",
       main = paste0(outGene_Amyg$Symbol[i], "\n", outGene_Amyg$gencodeID[i]) )
 points(gAmygExprs[i,] ~ jitter(as.numeric(factor(pdAmyg$Dx))),
       pch=21, bg=as.numeric(pdAmyg$Dx), cex=1.5)
 axis(1, at=1:1:length(levels(pdAmyg$Dx)), labels = levels(as.factor(pdAmyg$Dx)))
 legend("top", paste0("p=",signif(outGene_Amyg$P_PrimaryDxControl[i],3)), bg="white")
}
dev.off()


sigOrderMat = as.data.frame(apply(outExon_Amyg[,c(23:24)], 2, function(x) order(x)[1:100]))
ooL = sigOrderMat$"P_PrimaryDxControl"
pdf("plots/top_exons_MDD_vs_cnt_Amyg.pdf",h=6,w=6)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(4,"Dark2"))
for(i in ooL) {
  boxplot(eAmygExprs[i,] ~ as.numeric(pdAmyg$Dx) , outline=FALSE, xaxt="n",
	   pch = 21,
       cex=2, xlab="",
	   ylim = c(min(eAmygExprs[i,]),max(eAmygExprs[i,])),
       ylab="Residualized Expression",
       main = paste0(outExon_Amyg$Symbol[i], "\n", outExon_Amyg$coord[i]) )
 points(eAmygExprs[i,] ~ jitter(as.numeric(factor(pdAmyg$Dx))),
       pch=21, bg=as.numeric(pdAmyg$Dx), cex=1.5)
 axis(1, at=1:1:length(levels(pdAmyg$Dx)), labels = levels(as.factor(pdAmyg$Dx)))
 legend("top", paste0("p=",signif(outExon_Amyg$P_PrimaryDxControl[i],3)), bg="white")
}
dev.off()


sigOrderMat = as.data.frame(apply(outJxn_Amyg_anno[,c(25:26)], 2, function(x) order(x)[1:100]))
ooL = sigOrderMat$"P_PrimaryDxControl"
pdf("plots/top_jxns_MDD_vs_cnt_Amyg.pdf",h=6,w=6)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(4,"Dark2"))
for(i in ooL) {
  jxnI = rownames(outJxn_Amyg_anno)[i]
  indI = which(rownames(jAmygExprs) == jxnI)
  boxplot(jAmygExprs[indI,] ~ as.numeric(pdAmyg$Dx) , outline=FALSE, xaxt="n",
	   pch = 21,
       cex=2, xlab="",
	   ylim = c(min(jAmygExprs[indI,]),max(jAmygExprs[indI,])),
       ylab="Residualized Expression",
       main = paste0(outJxn_Amyg_anno$newGeneSymbol[i],", ",outJxn_Amyg_anno$Class[i],"\n", rownames(outJxn_Amyg_anno)[i]) )
 points(jAmygExprs[indI,] ~ jitter(as.numeric(factor(pdAmyg$Dx))),
       pch=21, bg=as.numeric(pdAmyg$Dx), cex=1.5)
 axis(1, at=1:1:length(levels(pdAmyg$Dx)), labels = levels(as.factor(pdAmyg$Dx)))
 legend("top", paste0("p=",signif(outJxn_Amyg_anno$P_PrimaryDxControl[i],3)), bg="white")
}
dev.off()


sigOrderMat = as.data.frame(apply(outTx_Amyg[,c(11:12)], 2, function(x) order(x)[1:100]))
ooL = sigOrderMat$"P_PrimaryDxControl"
pdf("plots/top_tx_MDD_vs_cnt_Amyg.pdf",h=6,w=6)
par(mar=c(5,6,5,2),cex.axis=2,cex.lab=2,cex.main=2)
palette(brewer.pal(4,"Dark2"))
for(i in ooL) {
  boxplot(tAmygExprs[i,] ~ as.numeric(pdAmyg$Dx) , outline=FALSE, xaxt="n",
	   pch = 21,
       cex=2, xlab="",
	   ylim = c(min(tAmygExprs[i,]),max(tAmygExprs[i,])),
       ylab="Residualized Expression",
       main = paste0(outTx_Amyg$gene_name[i], "\n", rownames(outTx_Amyg)[i]) )
 points(tAmygExprs[i,] ~ jitter(as.numeric(factor(pdAmyg$Dx))),
       pch=21, bg=as.numeric(pdAmyg$Dx), cex=1.5)
 axis(1, at=1:1:length(levels(pdAmyg$Dx)), labels = levels(as.factor(pdAmyg$Dx)))
 legend("top", paste0("p=",signif(outTx_Amyg$P_PrimaryDxControl[i],3)), bg="white")
}
dev.off()










## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
