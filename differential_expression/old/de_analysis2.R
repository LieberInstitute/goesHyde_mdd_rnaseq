#####################################
#Make qSVs (adapted from https://github.com/LieberInstitute/qsva_brain/blob/master/brainseq_phase2_qsv/make_qSVs.R)
#####################################

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library('readxl')
library('devtools')

setwd('/dcl02/lieber/ptsd/RNAseq/qsva/dlpfc_amygdala/')

#load rse_gene object
load('/dcl02/lieber/ptsd/RNAseq/count_data/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')

colData(rse_gene)$RIN = sapply(colData(rse_gene)$RIN,"[",1)
colData(rse_gene)$totalAssignedGene = sapply(colData(rse_gene)$totalAssignedGene, mean)
colData(rse_gene)$mitoRate = sapply(colData(rse_gene)$mitoRate,mean)
colData(rse_gene)$overallMapRate = sapply(colData(rse_gene)$overallMapRate, mean)
colData(rse_gene)$rRNA_rate = sapply(colData(rse_gene)$rRNA_rate,mean)
colData(rse_gene)$ERCCsumLogErr = sapply(colData(rse_gene)$ERCCsumLogErr,mean)


## model
with(colData(rse_gene), table(Group))

#> with(colData(rse_gene), table(Group))
#Group
#Control     MDD    PTSD 
#    431     432     422 


mod = model.matrix(~Group + AgeDeath + Sex + Race + Region + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate, data = colData(rse_gene))

#load cov_rse object from brainseq data
load("/dcl02/lieber/ptsd/RNAseq/qsva/dlpfc_amygdala/rdas/degradation_rse_PTSD_usingJoint.rda", verbose = TRUE)

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod)
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]

#> getPcaVars(qsvBonf)[1:k]
# [1] 75.500  5.300  3.170  1.760  1.230  1.020  0.816  0.407  0.377  0.352
#[11]  0.326  0.290  0.269  0.242  0.217  0.207  0.189  0.161  0.154

modQsva = cbind(mod, qSVs)

save(qsvBonf, qSVs, mod, modQsva, file = 'rdas/PTSD_qsvs.Rdata')
dir.create("pdf", showWarnings = FALSE)
pdf(file='pdf/qsvs_var_explained.pdf', useDingbats = FALSE)
plot(getPcaVars(qsvBonf)[1:k], pch=20)
dev.off()



#####################################
#Explore qSVs (adapted from https://github.com/LieberInstitute/qsva_brain/blob/master/brainseq_phase2_qsv/explore_qsvs.R)
#####################################

library(jaffelab)
library(rtracklayer)
library(recount)
library(recount.bwtool)
library(BiocParallel)
library(SummarizedExperiment)
library('readxl')
library('devtools')

setwd('/dcl02/lieber/ptsd/RNAseq/qsva/dlpfc_amygdala/')

# gene-level expression for degradation data
load('/dcl02/lieber/ptsd/RNAseq/count_data/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata')
rse_gene <- rse_gene[, rse_gene$Region %in% c('DLPFC', 'BasoAmyg', 'MedialAmyg', 'dACC')]

getRPKM = recount::getRPKM
geneRpkm = getRPKM(rse_gene, length="Length")
pd = colData(rse_gene)

# do pca on genes
pca = prcomp(t(log2(geneRpkm+1)))

# percent of variance explained by pcas
pcaVars = getPcaVars(pca)
getPcaVars(pca)[1:5]
# > getPcaVars(pca)[1:5]
#[1] 19.50 11.20  6.49  4.93  2.59


## degradation plots
pdf('pdf/degradation.pdf', useDingbats = FALSE)
plot(pca$x[,1] ~ pd$totalAssignedGene,
     xlab = "Gene Assignment Rate", pch=20,
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(pca$x[,1] ~ factor(pd$Region),
     xlab = "Region",
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"))
plot(pca$x[,2] ~ factor(pd$Region),
     xlab = "Region",
     ylab=paste0("pca2: ",getPcaVars(pca)[2],"% Var Expl"))
plot(pca$x[,1] ~ pd$RIN,
     xlab = "RIN", pch=20,
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(pca$x[,1] ~ pd$mitoRate,
     xlab = "mitoRate", pch=20,
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(pca$x[,1] ~ pd$rRNA_rate,
     xlab = "rRNA Rate", pch=20,
   ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(pca$x[,1] ~ pd$overallMapRate,
     xlab = "Overall Map Rate", pch=20,
   ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
dev.off()



#PTSD plots
pd$Group = factor(pd$Group)
pd$Sex = factor(pd$Sex)
pd$Race = factor(pd$Race)

pdf('pdf/PTSD_plots.pdf', useDingbats = FALSE)

plot(pca$x[,1] ~ pd$Group,
     xlab = "Group",
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"))

plot(pca$x[,1] ~ pd$Sex,
     xlab = "Sex",
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"))

plot(pca$x[,1] ~ pd$AgeDeath,
     xlab = "Age", pch=20,
     ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')

plot(pca$x[,2] ~ pd$AgeDeath,
     xlab = "Age", pch=20,
     ylab=paste0("pca2: ",getPcaVars(pca)[2],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')

plot(pca$x[,1] ~ pd$Race,
     xlab = "Race", pch=20,
	 ylab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"))


plot(pca$x[,2] ~ pd$Race,
     xlab = "Race", pch=20,
	 ylab=paste0("pca2: ",getPcaVars(pca)[2],"% Var Expl"))


dev.off()


# top 1000 expressed regions associated with degradation
load("rdas/degradation_rse_PTSD_usingJoint.rda", verbose = TRUE)

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))
getPcaVars(qsvBonf)[1:5]
# > getPcaVars(qsvBonf)[1:5]
#[1] 75.50  5.30  3.17  1.76  1.23

##PTSD qsv plots
pdf("pdf/qSVs_PTSD.pdf", useDingbats = FALSE)
plot(qsvBonf$x[,1] ~ pd$totalAssignedGene,
     xlab = "Gene Assignment Rate",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(qsvBonf$x[,1] ~ pd$Group,
     xlab = "Group",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
plot(qsvBonf$x[,1] ~ factor(pd$Region),
     xlab = "Region",	ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"))
plot(qsvBonf$x[,1] ~ pd$RIN,
     xlab = "RIN",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(qsvBonf$x[,1] ~ pd$mitoRate,
     xlab = "mitoRate",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(qsvBonf$x[,2] ~ pd$mitoRate,
     xlab = "mitoRate",pch=20,
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(qsvBonf$x[,1] ~ pd$AgeDeath,
     xlab = "Age",pch=20,
     ylab=paste0("qSV1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(qsvBonf$x[,2] ~ pd$AgeDeath,
     xlab = "Age",pch=20,
     ylab=paste0("qSV2: ",getPcaVars(qsvBonf)[2],"% Var Expl"), col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
legend('bottomright', sort(unique(pd$Region)), lwd=2, col = c("orange", "cadetblue2","magenta","black"), bty='n')
dev.off()


print('Regression: qsv1 vs Group + AgeDeath + Sex + Race + Region + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate')
summary(lm(qsvBonf$x[,1] ~ pd$Group + pd$AgeDeath + pd$Sex + pd$Race + pd$Region + pd$mitoRate + pd$rRNA_rate + pd$totalAssignedGene + pd$RIN + pd$overallMapRate))
#> print('Regression: qsv1 vs Group + AgeDeath + Sex + Race + Region + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate')
#[1] "Regression: qsv1 vs Group + AgeDeath + Sex + Race + Region + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate"
#> summary(lm(qsvBonf$x[,1] ~ pd$Group + pd$AgeDeath + pd$Sex + pd$Race + pd$Region + pd$mitoRate + pd$rRNA_rate + pd$totalAssignedGene + pd$RIN + pd$overallMapRate))
#
#Call:
#lm(formula = qsvBonf$x[, 1] ~ pd$Group + pd$AgeDeath + pd$Sex + 
#    pd$Race + pd$Region + pd$mitoRate + pd$rRNA_rate + pd$totalAssignedGene + 
#    pd$RIN + pd$overallMapRate)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-37.423  -4.598   0.716   5.308  24.289 
#
#Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          -9.614e+01  4.155e+00 -23.140  < 2e-16 ***
#pd$GroupMDD          -1.095e+00  5.569e-01  -1.966  0.04957 *  
#pd$GroupPTSD         -1.792e+00  5.882e-01  -3.047  0.00236 ** 
#pd$AgeDeath          -3.819e-02  1.663e-02  -2.296  0.02183 *  
#pd$SexM               1.534e+00  4.768e-01   3.217  0.00133 ** 
#pd$RaceCAUC          -8.534e-01  5.782e-01  -1.476  0.14018    
#pd$RaceHISP           1.966e+00  1.803e+00   1.091  0.27563    
#pd$RegiondACC         3.152e+00  6.516e-01   4.838 1.47e-06 ***
#pd$RegionDLPFC        2.600e+00  6.972e-01   3.730  0.00020 ***
#pd$RegionMedialAmyg  -8.334e-01  6.237e-01  -1.336  0.18172    
#pd$mitoRate           5.206e+01  2.220e+01   2.345  0.01920 *  
#pd$rRNA_rate         -4.612e+03  1.537e+03  -3.001  0.00275 ** 
#pd$totalAssignedGene  2.051e+02  3.941e+00  52.027  < 2e-16 ***
#pd$RIN                8.708e+00  2.916e-01  29.863  < 2e-16 ***
#pd$overallMapRate    -7.114e+01  4.117e+00 -17.278  < 2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 7.853 on 1270 degrees of freedom
#Multiple R-squared:  0.8478,	Adjusted R-squared:  0.8462 
#F-statistic: 505.5 on 14 and 1270 DF,  p-value: < 2.2e-16


#####################################
#Case-Control DLPFC (https://github.com/LieberInstitute/qsva_brain/blob/master/brainseq_phase2_qsv/casectrl_DLPFC.R)
#####################################

library(jaffelab)
library(SummarizedExperiment)
library(limma)
library(edgeR)
library('devtools')

setwd('/dcl02/lieber/ptsd/RNAseq/qsva/dlpfc_amygdala/')

#load rse object from PTSD data
load('/dcl02/lieber/ptsd/RNAseq/count_data/rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata', verbose = TRUE)


#load qsvBonf, qSVs, mod, modQsva object from PTSD data
load("/dcl02/lieber/ptsd/RNAseq/qsva/dlpfc_amygdala/rdas/PTSD_qsvs.Rdata", verbose = TRUE)

#filter for DLPFC
keepIndex = which(rse_gene$Region == "DLPFC")
rse_gene <- rse_gene[, keepIndex]

#Get rid of region columns
deleteIndex <- c('RegionDLPFC','RegiondACC','RegionMedialAmyg')
mod <- mod[keepIndex, !colnames(mod) %in% deleteIndex, drop = FALSE]
modQsva <- modQsva[keepIndex, !colnames(modQsva) %in% deleteIndex, drop = FALSE]

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
#calculate library-size adjustment
dge = calcNormFactors(dge)

pdf('pdf/dlpfc_voom_qsva.pdf', useDingbats = FALSE)
vGene = voom(dge,modQsva, plot=TRUE)
dev.off()
fitGene = lmFit(vGene)
eBGene = eBayes(fitGene)
sigGene = topTable(eBGene,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene = sigGene[rownames(rse_gene),]

## no qSVA
pdf('pdf/dlpfc_voom_noqsva.pdf', useDingbats = FALSE)
vGene0 = voom(dge,mod, plot=TRUE)
dev.off()
fitGene0 = lmFit(vGene0)
eBGene0 = eBayes(fitGene0)
sigGene0 = topTable(eBGene0,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGene0 = sigGene0[rownames(rse_gene),]

## no adjustment vars
pdf('pdf/dlpfc_voom_noadj.pdf', useDingbats = FALSE)
vGeneNoAdj = voom(dge, with(colData(rse_gene), model.matrix( ~ Group)), plot=TRUE)
dev.off()
fitGeneNoAdj = lmFit(vGeneNoAdj)
eBGeneNoAdj = eBayes(fitGeneNoAdj)
sigGeneNoAdj = topTable(eBGeneNoAdj,coef=2,
	p.value = 1,number=nrow(rse_gene))
outGeneNoAdj = sigGeneNoAdj[rownames(rse_gene),]

stopifnot(identical(rownames(outGene), rownames(outGene0)))
stopifnot(identical(rownames(outGene), rownames(outGeneNoAdj)))

save(outGene, outGene0, outGeneNoAdj, file = "rdas/groupStats_dlpfc_filtered_qSVA_geneLevel.rda")

write.csv(outGene, file="outGene.csv")
write.csv(outGene0, file="outGene0.csv")
write.csv(outGeneNoAdj, file="outGeneNoAdj.csv")



#####################################
#Explore Case-Control DLPFC (https://github.com/LieberInstitute/qsva_brain/blob/master/brainseq_phase2_qsv/explore_case_control.R)
#####################################







## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
