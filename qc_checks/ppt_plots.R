##

library(jaffelab)
library(SummarizedExperiment)
library(recount)
library(readxl)
library(RColorBrewer)
library(sessioninfo)
library(dplyr)
## load phenotype and alignment data
load("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/data/rse_gene_raw_GoesZandi_n1140.rda", verbose = TRUE) 
rse_gene <- rse_gene[,rse_gene$Experiment == "psychENCODE_MDD"]
pd_mdd <- colData(rse_gene) %>% as.data.frame

pd <- read.csv("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/read_and_alignment_metrics_goesHyde_MDD.csv", stringsAsFactors=FALSE, row.names=1)
pd <- pd[,!colnames(pd) %in% colnames(pd_mdd)]
pd$RNum <- rownames(pd)
pd <- left_join(pd_mdd, pd,by ='RNum')
rownames(pd) <- pd$RNum

table(pd$BrainRegion, pd$PrimaryDx)
table(pd$Sex, pd$PrimaryDx)
table(pd$Race, pd$PrimaryDx)
table(table(pd$BrNum))
summary(pd$Age)

# > table(pd$BrainRegion, pd$PrimaryDx)
# 
# Control MDD
# Amygdala      73 240
# sACC          71 239
# > table(pd$BrainRegion, pd$PrimaryDx)
# 
# Control MDD
# Amygdala      73 240
# sACC          71 239
# > table(pd$Sex, pd$PrimaryDx)
# 
# Control MDD
# F      35 159
# M     109 320
# > table(pd$Race, pd$PrimaryDx)
# 
# Control MDD
# AA         4   0
# AS         2   0
# CAUC     134 479
# HISP       4   0
# > table(table(pd$BrNum))
# 
# 1   2   3 
# 16 296   5 
# > summary(pd$Age)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 17.94   32.10   47.42   46.09   55.70   95.27


summary(pd$numReads)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 10284026  87697208 103523346 122609239 138632467 453612210 


aL=min(pd$totalAssignedGene)
aH=max(pd$totalAssignedGene)
mL=min(pd$overallMapRate)
mH=max(pd$overallMapRate)
mitoL=min(pd$rRNA_rate)
mitoH=max(pd$rRNA_rate)
rL=min(pd$numReads)
rH=max(pd$numReads)

# drop samples?
pdf("RIN_check_predrop.pdf", h=10,w=10)
par(mfcol=c(2,2),mar=c(5,6,2,2),cex.axis=1.8,cex.lab=1.8)
plot(pd$RIN, pd$overallMapRate, pch =21, bg="grey")
abline(h=0.5, lty=2)
legend("bottomright", "a)", bty="n", cex=2)
plot(pd$RIN, pd$totalAssignedGene, pch =21, bg="grey")
abline(h=0.3, lty=2)
legend("bottomright", "b)", bty="n", cex=2)
plot(pd$rRNA_rate, pd$overallMapRate, pch =21, bg="grey")
abline(h=0.5, lty=2)
abline(v=5e-4, lty=2)
legend("bottomright", "c)", bty="n", cex=2)
plot(pd$RIN, log10(pd$numReads), pch =21, bg="grey")
abline(h=log10(9e6), lty=2)
legend("bottomright", "d)", bty="n", cex=2)
dev.off()


####################################

## drop samples
pd$dropMetrics = FALSE
pd$dropMetrics[pd$overallMapRate < 0.5 | pd$totalAssignedGene < .3 | pd$numReads < 1e7] = TRUE

table(pd$dropMetrics)
# FALSE  TRUE 
# 618     5

####################################

pdf("RIN_check_postdrop.pdf", h=10,w=10)
par(mfcol=c(2,2),mar=c(5,6,2,2),cex.axis=1.8,cex.lab=1.8)
plot(pd$RIN, pd$overallMapRate, pch=21, bg=pd$dropMetrics+1, cex=pd$dropMetrics+1, ylim=c(mL,mH) )
abline(h=0.5, lty=2)
plot(pd$RIN, pd$totalAssignedGene, pch=21, bg=pd$dropMetrics+1, cex=pd$dropMetrics+1, ylim=c(aL,aH))
abline(h=0.3, lty=2)
plot(pd$rRNA_rate, pd$overallMapRate, pch=21, bg=pd$dropMetrics+1, cex=pd$dropMetrics+1, 
     xlim=c(mitoL,mitoH), ylim=c(mL,mH))
abline(h=0.5, lty=2)
plot(pd$RIN, log10(pd$numReads), pch=21, bg=pd$dropMetrics+1, cex=pd$dropMetrics+1, 
     ylim=c(log10(rL),log10(rH)))
abline(h=log10(1e7), lty=2)
dev.off()

##################################################

## actually drop the 3 low read samples, these have problems later
## none below 1e7 on n1140
pd = pd[pd$numReads > 1e7,]

##################################################

## metrics by flowcell

table(pd$Plate)	# Plate 2 are the 96 rerun samples

### plot
pdf("metrics_by_plate.pdf", h=5,w=5)
par(mar=c(7,6,2,2),cex.axis=1,cex.lab=1.5,cex.main=2)
palette(brewer.pal(8,"Dark2"))
boxplot(pd$overallMapRate ~ pd$Plate, las= 3,
        ylim = c(0,1), xlab="Plate",
        outline=FALSE,ylab="Overall Map Rate")
points(pd$overallMapRate ~ jitter(as.numeric(
  factor(pd$Plate)),
  amount=0.15), pch=21,
  bg = as.numeric(factor(pd$BrainRegion)))
legend("bottomleft", levels(factor(pd$BrainRegion)),
       pch = 15, col = 1:2,cex=.8)
abline(h=.5, lty=2, col="grey")

boxplot(pd$totalAssignedGene ~ pd$Plate, las= 3,
        ylim = c(0,.7), xlab="Plate",
        outline=FALSE,ylab="Gene Assignement Rate")
points(pd$totalAssignedGene ~ jitter(as.numeric(
  factor(pd$Plate)),
  amount=0.15), pch=21,
  bg = as.numeric(factor(pd$BrainRegion)))
abline(h=.3, lty=2, col="grey")

dev.off()


#################################
### regional labeling
library(rtracklayer)
gRpkm = recount::getRPKM(rse_gene, "Length")

## filter low expressing genes
gRpkm = gRpkm[which(rowMeans(gRpkm)>0.1),]
yExprs = log2(gRpkm+1)

### top 1000 genes different between regions
library(genefilter)
ngenes = 100

p = rowttests(yExprs, as.factor(pd$BrainRegion))
pOrd = p[order(p$p.value)[1:ngenes],]
ind1000 = which(rownames(p) %in% rownames(pOrd))
yExprs1000 = yExprs[ind1000,]

### estimate brain region
amyg = ifelse(pd$BrainRegion=="Amygdala",1,0)
sacc = ifelse(pd$BrainRegion=="sACC",1,0)
mod = data.frame(model.matrix(~amyg+sacc - 1))

library(limma)
fit = lmFit(yExprs1000, mod)
Xmat = fit$coef
Dmat = t(Xmat)%*%Xmat
guess = apply(yExprs1000, 2, function(x)  solve(Dmat, t(Xmat) %*% x))[2,]
stopifnot(identical(names(guess),paste0(rownames(pd),"_psychENCODE_MDD")))
pd$guess = guess

### plot
pdf("region_check_100.pdf", h=6,w=6)
par(mar=c(8,6,2,2),cex.axis=1.5,cex.lab=1.5,cex.main=2)
palette(brewer.pal(8,"Dark2"))
boxplot(pd$guess ~ pd$BrainRegion, las= 3,
        ylim = range(pd$guess),
        outline=FALSE,ylab="sACC Identity", xlab="")
points(pd$guess ~ jitter(as.numeric(
  factor(pd$BrainRegion)),
  amount=0.25), pch=21,
  bg = as.numeric(as.factor(pd$BrainRegion)))
segments(0,0.9,1.5,0.9, lty=2, col="grey")
segments(1.5,0.3,2.5,0.3, lty=2, col="grey")
dev.off()

info_cols <- c("SAMPLE_ID","BrNum","Age","Sex","Race","PrimaryDx","BrainRegion","dropMetrics","guess") #previoulsy included "dropGeno"
pd[which(pd$BrainRegion=="sACC" & pd$guess<.3),info_cols]
pd[which(pd$BrainRegion=="Amygdala" & pd$guess>0.9),info_cols]


# > pd[which(pd$BrainRegion=="sACC" & pd$guess<.3),info_cols]
# SAMPLE_ID  BrNum      Age Sex Race PrimaryDx BrainRegion dropMetrics        guess
# 9      R14179 Br1469 28.57000   M CAUC   Control        sACC       FALSE -0.008304666
# 49     R17520 Br1635 52.31000   M CAUC       MDD        sACC       FALSE  0.218813100
# 55     R17527 Br1675 32.10000   M CAUC       MDD        sACC       FALSE -0.077176996
# 74     R17547 Br1754 25.98000   F CAUC       MDD        sACC       FALSE  0.012868825
# 106    R17579 Br6099 26.66381   M CAUC       MDD        sACC       FALSE  0.155488347
# 407    R17894 Br8017 49.29500   M CAUC       MDD        sACC       FALSE  0.274937417
# 454    R17941 Br8133 62.36813   M CAUC       MDD        sACC       FALSE -0.132722760
# 464    R17951 Br2333 68.85000   F CAUC   Control        sACC       FALSE -0.051220201
# 465    R17952 Br5549 67.47000   M CAUC       MDD        sACC       FALSE  0.012867425
# 489    R18430 Br3863 69.33904   M CAUC   Control        sACC       FALSE -0.133637833
# 514    R18458 Br8313 55.41672   M CAUC   Control        sACC       FALSE  0.140951403
# 619    R19186 Br1475 43.30000   M CAUC       MDD        sACC       FALSE  0.081657998
# > pd[which(pd$BrainRegion=="Amygdala" & pd$guess>0.9),info_cols]
# SAMPLE_ID  BrNum   Age Sex Race PrimaryDx BrainRegion dropMetrics    guess
# 26    R17496 Br1469 28.57   M CAUC   Control    Amygdala       FALSE 1.073439

## Br1469 labels got switched
pd["R14179","BrainRegion"] = "Amygdala"
pd["R17496","BrainRegion"] = "sACC"

## Drop the rest
dropInd = which( (pd$BrainRegion=="sACC" & pd$guess<.3) |
                   (pd$BrainRegion=="Amygdala" & pd$guess>0.9)  )
pd$dropRegion = FALSE
pd$dropRegion[dropInd] = TRUE

table(pd$dropRegion)
# FALSE  TRUE
# 612    11

########################################################.
gia <- read.csv("../synapse/genotypeInferredAncestry.csv")
pd$geneticRace = pd$Race
pd$geneticRace[match(gia$BrNum, pd$BrNum)] = gia$genotypeInferredAncestry

pd$dropRace = pd$geneticRace != "CAUC"

table(pd$Race, pd$geneticRace)
#       AA ambiguous  AS CAUC HISP
# AA     4         0   0    0    0
# AS     0         0   2    0    0
# CAUC   0         2   0  611    0
# HISP   0         0   0    1    3

table(pd$dropRace)
# FALSE  TRUE 
# 612    11 


########################################################.
info_cols <- c("SAMPLE_ID","BrNum","Age","Sex","Race","PrimaryDx","BrainRegion","RIN","Plate","dropMetrics","dropRegion","dropRace") #previoulsy included "dropGeno"
pd[which(pd$dropMetrics==TRUE | pd$dropRegion==TRUE | pd$dropRace==TRUE),info_cols]

# SAMPLE_ID  BrNum      Age Sex Race PrimaryDx BrainRegion RIN Plate overallMapRate totalAssignedGene dropMetrics dropRegion dropRace
# 9      R14179 Br1469 28.57000   M CAUC   Control        sACC 7.6     2         0.8547         0.4880644       FALSE       TRUE    FALSE
# 26     R17496 Br1469 28.57000   M CAUC   Control    Amygdala 8.0     2         0.8149         0.4895696       FALSE       TRUE    FALSE
# 49     R17520 Br1635 52.31000   M CAUC       MDD        sACC 6.6     6         0.8219         0.4278223       FALSE       TRUE    FALSE
# 55     R17527 Br1675 32.10000   M CAUC       MDD        sACC 6.8     2         0.7848         0.4126576       FALSE       TRUE    FALSE
# 66     R17538 Br1723 55.70000   F CAUC       MDD    Amygdala 5.3     1         0.7946         0.2045859        TRUE      FALSE    FALSE
# 74     R17547 Br1754 25.98000   F CAUC       MDD        sACC 5.7     5         0.7315         0.4009698       FALSE       TRUE    FALSE
# 106    R17579 Br6099 26.66381   M CAUC       MDD        sACC 6.9     7         0.8584         0.4527738       FALSE       TRUE    FALSE
# 140    R17616 Br1986 25.35000   M CAUC       MDD    Amygdala 6.7     5         0.8204         0.4589626       FALSE      FALSE     TRUE
# 311    R17794 Br5763 53.89000   M CAUC       MDD        sACC 7.1     1         0.4800         0.4707335        TRUE      FALSE    FALSE
# 354    R17839 Br5965 30.12183   M CAUC       MDD    Amygdala 5.6     4         0.7764         0.2037251        TRUE      FALSE    FALSE
# 393    R17880 Br6227 18.62822   M HISP   Control        sACC 7.4     5         0.9070         0.4828394       FALSE      FALSE     TRUE
# 407    R17894 Br8017 49.29500   M CAUC       MDD        sACC 7.2     3         0.7022         0.3977602       FALSE       TRUE    FALSE
# 427    R17914 Br8053 47.00616   M CAUC   Control    Amygdala 5.5     1         0.2259         0.4183653        TRUE      FALSE    FALSE
# 454    R17941 Br8133 62.36813   M CAUC       MDD        sACC 5.7     1         0.7399         0.4075764       FALSE       TRUE    FALSE
# 464    R17951 Br2333 68.85000   F CAUC   Control        sACC 5.0     4         0.7008         0.4349815       FALSE       TRUE    FALSE
# 465    R17952 Br5549 67.47000   M CAUC       MDD        sACC 5.6     5         0.7961         0.4222974       FALSE       TRUE    FALSE
# 476    R17963 Br5526 68.85000   M CAUC       MDD    Amygdala 5.4     4         0.1883         0.4040689        TRUE      FALSE    FALSE
# 487    R18425 Br6285 64.82000   M   AS   Control    Amygdala 6.8     3         0.8767         0.4339124       FALSE      FALSE     TRUE
# 488    R18426 Br6285 64.82000   M   AS   Control        sACC 7.2     3         0.8286         0.4597279       FALSE      FALSE     TRUE
# 489    R18430 Br3863 69.33904   M CAUC   Control        sACC 6.5     7         0.7480         0.4311182       FALSE       TRUE    FALSE
# 490    R18431 Br3872 59.24162   F HISP   Control    Amygdala 5.8     4         0.6814         0.4078292       FALSE      FALSE     TRUE
# 491    R18432 Br3872 59.24162   F HISP   Control        sACC 5.9     4         0.7363         0.3812791       FALSE      FALSE     TRUE
# 514    R18458 Br8313 55.41672   M CAUC   Control        sACC 5.5     4         0.8069         0.4208956       FALSE       TRUE    FALSE
# 515    R18459 Br6123 30.76249   F   AA   Control    Amygdala 6.6     3         0.8757         0.4761234       FALSE      FALSE     TRUE
# 516    R18460 Br6123 30.76249   F   AA   Control        sACC 7.1     4         0.8869         0.4457243       FALSE      FALSE     TRUE
# 517    R18461 Br5146 47.01000   F   AA   Control    Amygdala 6.8     4         0.8658         0.4538249       FALSE      FALSE     TRUE
# 518    R18462 Br5146 47.01000   F   AA   Control        sACC 7.6     4         0.9031         0.4928282       FALSE      FALSE     TRUE
# 531    R18824 Br3877 51.85775   M CAUC       MDD    Amygdala 6.4     6         0.8369         0.4241977       FALSE      FALSE     TRUE
# 619    R19186 Br1475 43.30000   M CAUC       MDD        sACC 8.9     4         0.8739         0.4635724       FALSE       TRUE    FALSE

## plus 3 with really low reads that were also dropMetrics

pd$dropSum = rowSums(pd[,c("dropMetrics","dropRegion","dropRace")])
sum(pd$dropSum > 0) #27
# sum(pd$dropSum > 0) + 3
# # 41

qcresults = pd[,info_cols]
write.csv(qcresults, file="qc_dropping_results.csv")

pd = pd[-which(pd$dropSum > 0),]
nrow(pd)
# 596

table(pd$BrainRegion, pd$PrimaryDx)
table(pd$Sex, pd$PrimaryDx)
table(pd$Race, pd$PrimaryDx)
table(table(pd$BrNum))
summary(pd$Age)

# > table(pd$BrainRegion, pd$PrimaryDx)

# Control MDD
# Amygdala      68 235
# sACC          63 230
# > table(pd$Sex, pd$PrimaryDx)
# 
# Control MDD
# F      28 157
# M     103 308
# > table(pd$Race, pd$PrimaryDx)
# 
# Control MDD
# CAUC     130 465
# HISP       1   0
# > table(table(pd$BrNum))
# 
# 1   2   3 
# 31 275   5 
# > summary(pd$Age)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 17.94   32.17   47.33   45.98   54.65   95.27


message("Remaining Samples: n", nrow(qcresults))

##############
## PCA
rownames(pd) <- paste0(rownames(pd),"_psychENCODE_MDD")
## Filter for good sampels
gRpkm = recount::getRPKM(rse_gene, "Length")
gRpkm = gRpkm[,rownames(pd)]

## filter low expressing genes
gRpkm = gRpkm[which(rowMeans(gRpkm)>0.2),]
yExprs = log2(gRpkm+1)

pca1 = prcomp(t(yExprs))
pcaVars1 = getPcaVars(pca1)

pdf("pca_log2Rpkm_PC1_2.pdf", h=6,w=6)
par(mar=c(5,6,4,2),cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
palette(brewer.pal(8,"Spectral"))
plot(pca1$x, pch = 21, bg=as.numeric(factor(pd$PrimaryDx))*3,cex=2, main="Gene PCs",
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("topleft", paste0(levels(factor(pd$PrimaryDx))),
       pch = 15, col = c(3,6),cex=.9)
plot(pca1$x, pch = 21, bg=as.numeric(factor(pd$Sex))*2,cex=2, main="Gene PCs",
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("topleft", paste0(levels(factor(pd$Sex))),
       pch = 15, col = c(2,4),cex=.9)
plot(pca1$x, pch = 21, bg=as.numeric(factor(pd$BrainRegion))*3,cex=2, main="Gene PCs",
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("topleft", paste0(levels(factor(pd$BrainRegion))),
       pch = 15, col = c(3,6),cex=.9)
plot(pca1$x, pch = 21, bg=as.numeric(factor(pd$Race))*2,cex=2, main="Gene PCs",
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("topleft", paste0(levels(factor(pd$Race))),
       pch = 15, col = c(2,4,6,8),cex=.9)
plot(pca1$x, pch = 21, bg=as.numeric(factor(pd$Plate)),cex=2, main="Gene PCs",
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("topleft", paste0("Plate ",levels(factor(pd$Plate))),
       pch = 15, col = 1:8,cex=.9)
dev.off()




## PCA - add in bipolar samples
load("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/data/rse_gene_raw_GoesZandi_n1140.rda", verbose = TRUE) 
gRpkmBP = recount::getRPKM(rse_gene, "Length")
## filter low expressing genes
gRpkmBP = gRpkmBP[which(rowMeans(gRpkmBP)>0.2),]
yExprsBP = log2(gRpkmBP+1)

###### combine
# intersection of genes
# common <- intersect(rownames(yExprs), rownames(yExprsBP))
# yExprsComb = cbind(yExprs[common,],yExprsBP[common,])
yExprsComb = yExprsBP
# combined phenodata
pd = colData(rse_gene)[,c("BrNum","BrainRegion","PrimaryDx","Experiment")]
pd$group <- paste0(pd$Experiment,"_",pd$PrimaryDx)
pca1 = prcomp(t(yExprsComb))
pcaVars1 = getPcaVars(pca1)

pdf("pca_log2Rpkm_PC1_2_combined_datasets.pdf", h=6,w=6)
par(mar=c(5,6,4,2),cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
palette(brewer.pal(8,"Spectral"))
plot(pca1$x, pch = 21, bg=(as.numeric(factor(pd$Experiment))*4)-1,cex=2, main="Gene PCs",
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("topleft", paste0(levels(factor(pd$Experiment))),
       pch = 15, col = c(3,7),cex=.9)
plot(pca1$x, pch = 21, bg=as.numeric(factor(pd$PrimaryDx))*2,cex=2, main="Gene PCs",
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("topleft", paste0(levels(factor(pd$PrimaryDx))),
       pch = 15, col = c(2,4,6),cex=.9)
plot(pca1$x, pch = 21, bg=as.numeric(factor(pd$BrainRegion))*3,cex=2, main="Gene PCs",
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("topleft", paste0(levels(factor(pd$BrainRegion))),
       pch = 15, col = c(3,6),cex=.9)
plot(pca1$x, pch = 21, bg=as.numeric(factor(pd$group)),cex=2, main="Gene PCs",
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("topleft", paste0(levels(factor(pd$group))),
       pch = 15, col = 1:4,cex=.9)
dev.off()

# sgejobs::job_single('ppt_plots', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript ppt_plots.R")


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()


