##

library(jaffelab)
library(SummarizedExperiment)
library(recount)
library(readxl)
library(RColorBrewer)

## load phenotype and alignment data
pd <- read.csv("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/read_and_alignment_metrics_goesHyde_MDD.csv", stringsAsFactors=FALSE, row.names=1)
pd_mdd <- read.csv("../data/GoesMDD_pd_n1140.csv")
pd <- subset(pd, SAMPLE_ID %in% pd_mdd$RNum)

table(pd$BrainRegion, pd$PrimaryDx)
table(pd$Sex, pd$PrimaryDx)
table(pd$Race, pd$PrimaryDx)
table(table(pd$BrNum))
summary(pd$Age)

# > table(pd$BrainRegion, pd$PrimaryDx)
# 
# Control MDD
# Amygdala      73 240
# sACC          70 240
# > table(pd$Sex, pd$PrimaryDx)
# 
# Control MDD
# F      34 160
# M     109 320
# > table(pd$Race, pd$PrimaryDx)
# 
# Control MDD
# AA         4   0
# AS         2   0
# CAUC     133 480
# HISP       4   0
# > table(table(pd$BrNum))
# 
# 1   2 
# 7 308 
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

pd$Plate	# Plate 2 are the 96 rerun samples

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
### genotype check
### taken care of in MDD swap
# library(VariantAnnotation)
# library(rtracklayer)
# library(pheatmap)
# 
# ### Load in snpMap
# snpMap = rtracklayer::import("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Genotyping/common_missense_SNVs_hg38.bed")
# 
# ##########################
# # read in merged VCF file
# # point to your samples
# mergedVcfFile = '../preprocessed_data/Genotypes/mergedVariants.vcf.gz'
# vcf = readVcf(mergedVcfFile, genome=seqinfo(snpMap))
# info(vcf)$RS = mcols(snpMap)$name[match(rowRanges(vcf),snpMap)]		# may get warning about RS, can ignore
# 
# vcf = vcf[,match(pd$bamFile,rownames(colData(vcf)))]	#  subset to 623
# colnames(vcf) = pd$SAMPLE_ID
# 
# ######################
# # subset to high-depth
# vcf = vcf[info(vcf)$DP > 5*ncol(vcf) & info(vcf)$DP < 80 *ncol(vcf) &
#           nchar(ref(vcf)) == 1 & elementNROWS(alt(vcf)) == 1 &
#           info(vcf)$VDB >0.1,]
# 
# ########################################
# # plot snp correlation of all samples
# snps = geno(vcf)$GT
# snps[snps == "."] = 0
# snps[snps == "0/1"] = 1
# snps[snps == "1/1"] = 2
# class(snps) = "numeric"
# snpCor = cor(snps, use="pairwise.complete.obs")
# 
# ######################
# ## how do you want samples to be labeled in the plot?
# rownames(snpCor) = colnames(snpCor) = pd$BrNum
# 
# snpCor2 = snpCor[order(rownames(snpCor)),order(colnames(snpCor))]
# s = snpCor2[1:622,1:622]
# s[which(is.na(s), arr.ind=TRUE)] = 0
# s = round(s,3)
# 
# brInd = split(1:622, colnames(s))
# brMats = lapply(brInd, function(x) s[x,x])
# brMats = brMats[sapply(brMats, length)>1]
# 
# badMatch = brMats[sapply(brMats, function(x) length(which(x < 0.65))>1 )]
# names(badMatch)
#  # [1] "Br1582" "Br1698" "Br1910" "Br2084" "Br5486" "Br5615" "Br5694" "Br5930" "Br5956" "Br8148"
#  
# ## do those brains match any other brain the dataset instead?
# badInd = which(colnames(snpCor) %in% names(badMatch))
# snpCorBad = snpCor[badInd,]
# otherMatches = colnames(snpCorBad)[which(snpCorBad > 0.65 & snpCorBad < 1, arr.ind=TRUE)[,2]]
# 
# toPlot = which(rownames(snpCor) %in% c(names(badMatch), otherMatches))
# snpCorPlot = snpCor[toPlot,toPlot]
# 
# ######################
# ## plot
# library(RColorBrewer)
# col.pal = brewer.pal(9,"Blues")
# 
# pal = c(brewer.pal(12,"Paired"),"gold3","gray95","gray45","gray20")
# 
# # Data frame with column annotations.
# mat_col <- data.frame(Brain = unique(rownames(snpCorPlot)))
# rownames(mat_col) = mat_col$Brain
# 
# # List with colors for each annotation.
# mat_colors <- list(Brain = pal)
# names(mat_colors$Brain) <- unique(mat_col$Brain)
# 
# pdf("genotype_heatmap_new.pdf",h=10,w=10)
# pheatmap(snpCorPlot, 
# 		cluster_rows=T, 
# 		cluster_cols=T,
# 		color=col.pal,
# 		annotation_row = mat_col,
# 		annotation_colors = mat_colors)
# dev.off()
# 
### drop:
# dropBoth = c("Br1910","Br5956","Br5807")
# dropAmyg = c("Br5615","Br5930")
# dropsACC = c("Br5694","Br2084","Br5486")
# pd$dropGeno = FALSE
# pd$dropGeno[pd$BrNum %in% dropBoth] = TRUE
# pd$dropGeno[pd$BrNum %in% dropAmyg & pd$BrainRegion=="Amygdala"] = TRUE
# pd$dropGeno[pd$BrNum %in% dropsACC & pd$BrainRegion=="sACC"] = TRUE

# table(pd$dropGeno)
# # FALSE  TRUE
#   # 621    10



#################################
### regional labeling
library(rtracklayer)

load("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/rse_gene_goesHyde_MDD_n634.Rdata", verbose = TRUE)
gRpkm = recount::getRPKM(rse_gene, "Length")
gRpkm = gRpkm[,rownames(pd)]

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
stopifnot(identical(names(guess),pd$SAMPLE_ID))
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
# R14179    R14179 Br1469 28.57000   M CAUC   Control        sACC       FALSE -0.008304666
# R17520    R17520 Br1635 52.31000   M CAUC       MDD        sACC       FALSE  0.218813100
# R17527    R17527 Br1675 32.10000   M CAUC       MDD        sACC       FALSE -0.077176996
# R17547    R17547 Br1754 25.98000   F CAUC       MDD        sACC       FALSE  0.012868825
# R17579    R17579 Br6099 26.66381   M CAUC       MDD        sACC       FALSE  0.155488347
# R17894    R17894 Br8017 49.29500   M CAUC       MDD        sACC       FALSE  0.274937417
# R17941    R17941 Br8133 62.36813   M CAUC       MDD        sACC       FALSE -0.132722760
# R17951    R17951 Br5526 68.85000   M CAUC       MDD        sACC       FALSE -0.051220201
# R17952    R17952 Br5549 67.47000   M CAUC       MDD        sACC       FALSE  0.012867425
# R18430    R18430 Br3863 69.33904   M CAUC   Control        sACC       FALSE -0.133637833
# R18458    R18458 Br8313 55.41672   M CAUC   Control        sACC       FALSE  0.140951403
# R19186    R19186 Br1475 43.30000   M CAUC       MDD        sACC       FALSE  0.081657998
#
# > pd[which(pd$BrainRegion=="Amygdala" & pd$guess>0.9),info_cols]
# SAMPLE_ID  BrNum   Age Sex Race PrimaryDx BrainRegion dropMetrics    guess
# R17496    R17496 Br1469 28.57   M CAUC   Control    Amygdala       FALSE 1.073439

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
  # 620    11
  
########################################################.
gia <- read.csv("../synapse/genotypeInferredAncestry.csv")
pd$geneticRace = pd$Race
pd$geneticRace[match(gia$BrNum, pd$BrNum)] = gia$genotypeInferredAncestry

pd$dropRace = pd$geneticRace != "CAUC"

table(pd$dropRace)
# FALSE  TRUE 
# 613    10 


########################################################.
info_cols <- c("SAMPLE_ID","BrNum","Age","Sex","Race","PrimaryDx","BrainRegion","RIN","Plate","dropMetrics","dropRegion","dropRace") #previoulsy included "dropGeno"
pd[which(pd$dropMetrics==TRUE | pd$dropGeno==TRUE | pd$dropRegion==TRUE | pd$dropRace==TRUE),info_cols]

#         SAMPLE_ID  BrNum      Age Sex Race PrimaryDx BrainRegion RIN Plate dropMetrics dropRegion dropRace
# R17520    R17520 Br1635 52.31000   M CAUC       MDD        sACC 6.6     6       FALSE       TRUE    FALSE
# R17527    R17527 Br1675 32.10000   M CAUC       MDD        sACC 6.8     2       FALSE       TRUE    FALSE
# R17538    R17538 Br1723 55.70000   F CAUC       MDD    Amygdala 5.3     1        TRUE      FALSE    FALSE
# R17547    R17547 Br1754 25.98000   F CAUC       MDD        sACC 5.7     5       FALSE       TRUE    FALSE
# R17579    R17579 Br6099 26.66381   M CAUC       MDD        sACC 6.9     7       FALSE       TRUE    FALSE
# R17602    R17602 Br1910 25.61000   F CAUC       MDD    Amygdala 7.2     3       FALSE      FALSE    FALSE
# R17603    R17603 Br1910 25.61000   F CAUC       MDD        sACC 5.7     3       FALSE      FALSE    FALSE
# R17623    R17623 Br2084 48.07000   M CAUC   Control        sACC 6.2     2       FALSE      FALSE    FALSE
# R17734    R17734 Br5486 82.26000   M CAUC   Control        sACC 6.8     4       FALSE      FALSE    FALSE
# R17738    R17738 Br5694 59.94000   M CAUC       MDD        sACC 7.2     4       FALSE      FALSE    FALSE
# R17783    R17783 Br5615 18.93000   M CAUC       MDD    Amygdala 6.5     6       FALSE      FALSE    FALSE
# R17794    R17794 Br5763 53.89000   M CAUC       MDD        sACC 7.1     1        TRUE      FALSE    FALSE
# R17831    R17831 Br5930 49.44022   M CAUC       MDD    Amygdala 5.8     1       FALSE      FALSE    FALSE
# R17833    R17833 Br5956 57.51141   F CAUC       MDD    Amygdala 5.3     5       FALSE      FALSE    FALSE
# R17834    R17834 Br5956 57.51141   F CAUC       MDD        sACC 7.3     5       FALSE      FALSE    FALSE
# R17839    R17839 Br5965 30.12183   M CAUC       MDD    Amygdala 5.6     4        TRUE      FALSE    FALSE
# R17879    R17879 Br6227 18.62822   M HISP   Control    Amygdala 7.3     5       FALSE      FALSE     TRUE
# R17880    R17880 Br6227 18.62822   M HISP   Control        sACC 7.4     5       FALSE      FALSE     TRUE
# R17894    R17894 Br8017 49.29500   M CAUC       MDD        sACC 7.2     3       FALSE       TRUE    FALSE
# R17914    R17914 Br8053 47.00616   M CAUC   Control    Amygdala 5.5     1        TRUE      FALSE    FALSE
# R17941    R17941 Br8133 62.36813   M CAUC       MDD        sACC 5.7     1       FALSE       TRUE    FALSE
# R17951    R17951 Br5526 68.85000   M CAUC       MDD        sACC 5.0     4       FALSE       TRUE    FALSE
# R17952    R17952 Br5549 67.47000   M CAUC       MDD        sACC 5.6     5       FALSE       TRUE    FALSE
# R17953    R17953 Br5807 36.93000   F CAUC       MDD        sACC 6.4     7       FALSE      FALSE    FALSE
# R17963    R17963 Br5526 68.85000   M CAUC       MDD    Amygdala 5.4     4        TRUE      FALSE    FALSE
# R18425    R18425 Br6285 64.82000   M   AS   Control    Amygdala 6.8     3       FALSE      FALSE     TRUE
# R18426    R18426 Br6285 64.82000   M   AS   Control        sACC 7.2     3       FALSE      FALSE     TRUE
# R18430    R18430 Br3863 69.33904   M CAUC   Control        sACC 6.5     7       FALSE       TRUE    FALSE
# R18431    R18431 Br3872 59.24162   F HISP   Control    Amygdala 5.8     4       FALSE      FALSE     TRUE
# R18432    R18432 Br3872 59.24162   F HISP   Control        sACC 5.9     4       FALSE      FALSE     TRUE
# R18458    R18458 Br8313 55.41672   M CAUC   Control        sACC 5.5     4       FALSE       TRUE    FALSE
# R18459    R18459 Br6123 30.76249   F   AA   Control    Amygdala 6.6     3       FALSE      FALSE     TRUE
# R18460    R18460 Br6123 30.76249   F   AA   Control        sACC 7.1     4       FALSE      FALSE     TRUE
# R18461    R18461 Br5146 47.01000   F   AA   Control    Amygdala 6.8     4       FALSE      FALSE     TRUE
# R18462    R18462 Br5146 47.01000   F   AA   Control        sACC 7.6     4       FALSE      FALSE     TRUE
# R19186    R19186 Br1475 43.30000   M CAUC       MDD        sACC 8.9     4       FALSE       TRUE    FALSE

## plus 3 with really low reads that were also dropMetrics

pd$dropSum = rowSums(pd[,c("dropMetrics","dropRegion","dropRace")])
sum(pd$dropSum > 0) #26
# sum(pd$dropSum > 0) + 3
# # 41

pd = pd[-which(pd$dropSum > 0),]
nrow(pd)
# 597

table(pd$BrainRegion, pd$PrimaryDx)
table(pd$Sex, pd$PrimaryDx)
table(pd$Race, pd$PrimaryDx)
table(table(pd$BrNum))
summary(pd$Age)

# table(pd$BrainRegion, pd$PrimaryDx)
# 
# Control MDD
# Amygdala      67 237
# sACC          63 230
# > table(pd$Sex, pd$PrimaryDx)
# 
# Control MDD
# F      28 158
# M     102 309
# > table(pd$Race, pd$PrimaryDx)
# 
# Control MDD
# CAUC     130 467
# > table(table(pd$BrNum))
# 
# 1   2 
# 21 288 
# > summary(pd$Age)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 17.94   32.19   47.33   46.00   54.62   95.27
  
  
info_cols <- c("SAMPLE_ID","BrNum","Age","Sex","Race","PrimaryDx","BrainRegion","RIN","Plate","overallMapRate","totalAssignedGene","dropMetrics","dropRegion","dropRace")  
qcresults = pd[,info_cols]
write.csv(qcresults, file="qc_dropping_results.csv")




##############
## PCA


## of 593 samples
load("../preprocessed_data/rse_gene_goesHyde_MDD_n634.Rdata", verbose=TRUE)
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

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rse_gene_zandiHyde_Bipolar_LIBD_n540.Rdata", verbose=TRUE)
gRpkmBP = recount::getRPKM(rse_gene, "Length")
## filter low expressing genes
gRpkmBP = gRpkmBP[which(rowMeans(gRpkmBP)>0.2),]
yExprsBP = log2(gRpkmBP+1)

###### combine
# intersection of genes
common <- intersect(rownames(yExprs), rownames(yExprsBP))
yExprsComb = cbind(yExprs[common,],yExprsBP[common,])

# combined phenodata
pdBP = colData(rse_gene)[,1:12]
names(pdBP)[5] = "BrainRegion"
pdBP$BrNum = as.character(pdBP$BrNum)
pdBP$BrainRegion = as.character(pdBP$BrainRegion)
pdBP$PrimaryDx = as.character(pdBP$PrimaryDx)

pd$Exp = "GoesMDD"
pdBP$Exp = "ZandiBPD"
pdComb = rbind(pd[,c("BrNum","BrainRegion","PrimaryDx","Exp")], pdBP[,c("BrNum","BrainRegion","PrimaryDx","Exp")])
pdComb$group = paste0(pdComb$Exp,"_",pdComb$PrimaryDx)

pca1 = prcomp(t(yExprsComb))
pcaVars1 = getPcaVars(pca1)

pdf("pca_log2Rpkm_PC1_2_combined_datasets.pdf", h=6,w=6)
par(mar=c(5,6,4,2),cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
palette(brewer.pal(8,"Spectral"))
plot(pca1$x, pch = 21, bg=(as.numeric(factor(pdComb$Exp))*4)-1,cex=2, main="Gene PCs",
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("topleft", paste0(levels(factor(pdComb$Exp))),
       pch = 15, col = c(3,7),cex=.9)
plot(pca1$x, pch = 21, bg=as.numeric(factor(pdComb$PrimaryDx))*2,cex=2, main="Gene PCs",
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("topleft", paste0(levels(factor(pdComb$PrimaryDx))),
       pch = 15, col = c(2,4),cex=.9)
plot(pca1$x, pch = 21, bg=as.numeric(factor(pdComb$BrainRegion))*3,cex=2, main="Gene PCs",
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("topleft", paste0(levels(factor(pdComb$BrainRegion))),
       pch = 15, col = c(3,6),cex=.9)
plot(pca1$x, pch = 21, bg=as.numeric(factor(pdComb$group)),cex=2, main="Gene PCs",
     xlab=paste0("PC1: ", pcaVars1[1], "% Var Expl"),
     ylab=paste0("PC2: ", pcaVars1[2], "% Var Expl"))
legend("topleft", paste0(levels(factor(pdComb$group))),
       pch = 15, col = 1:8,cex=.9)
dev.off()


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# sgejobs::job_single('ppt_plots', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript ppt_plots.R")

