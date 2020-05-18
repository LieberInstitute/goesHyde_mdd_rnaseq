#####################################
# from /dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/R_scripts/DEG_scripts/qSV_model_analysis_DLPFC_PTSD.R
#####################################

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)
library(edgeR)
library(sessioninfo)

#setwd('/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/differential_expression/')

###############################
##### Load, clean, combine ####
###############################

#load objects
load('../exprs_cutoff/rse_gene.Rdata', verbose=TRUE)
load('../exprs_cutoff/rse_exon.Rdata', verbose=TRUE)
load('../exprs_cutoff/rse_jxn.Rdata', verbose=TRUE)
load('../exprs_cutoff/rse_tx.Rdata', verbose=TRUE)


pd = colData(rse_gene)

table(pd$Experiment)
 # GoesMDD ZandiBPD
     # 593      500

table(pd$BrainRegion)
# Amygdala     sACC
     # 541      552

table(pd$BrainRegion, pd$Sex)
            # F   M
  # Amygdala 163 378
  # sACC     166 386


table(pd$Experiment, pd$PrimaryDx)
           # MDD Control Bipolar
  # GoesMDD  463     130       0
  # ZandiBPD   0     255     245

table(pd$BrainRegion,pd$PrimaryDx)
           # MDD Control Bipolar
  # Amygdala 235     186     120
  # sACC     228     199     125


summary(pd$AgeDeath)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 17.37   34.44   47.21   46.57   55.86   95.27


##############################
############ load degredation
# #
# load("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/data/degradation_rse_MDDseq_BothRegions.rda")
# cov_rsemdd = cov_rse
# load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/data/degradation_rse_BipSeq_BothRegions.rda")
# cov_rsebip = cov_rse

load("../data/degradation_rse_MDDseq_BiPSeq_BothRegions.Rdata", verbose = TRUE)


# # combine
# cov_rsemdd$AgeDeath = cov_rsemdd$Age
# cov_rsemdd$RNum = cov_rsemdd$SAMPLE_ID
# cov_rsemdd$BrNum = as.character(cov_rsemdd$BrNum)
# colKeep = c("RNum","BrNum","Sex","Race","AgeDeath","BrainRegion","RIN","PrimaryDx","overallMapRate","totalAssignedGene","mitoRate","rRNA_rate")
# colData(cov_rsemdd) = colData(cov_rsemdd)[,colKeep]
# colData(cov_rsebip) = colData(cov_rsebip)[,colKeep]
#
# cov_rseboth = cbind(cov_rsemdd, cov_rsebip)
# cov_rseboth = cov_rseboth[,match(rse_gene$RNum,cov_rseboth$RNum)]	# put in order / drop dropped samples
# stopifnot(identical(cov_rseboth$RNum, rse_gene$RNum))
#
# cov_rse = cov_rseboth



## add ancestry
load("../genotype_data/goesHyde_bipolarMdd_Genotypes_n588_mds.rda", verbose = TRUE)

## keep samples with genotypes (remember that rse_gene is summarized experiment, so BrNum is in colData)
rse_gene <- rse_gene[, rse_gene$BrNum %in% rownames(mds)]
rse_exon <- rse_exon[, rse_exon$BrNum %in% rownames(mds)]
rse_jxn <- rse_jxn[, rse_jxn$BrNum %in% rownames(mds)]
rse_tx <- rse_tx[, rse_tx$BrNum %in% rownames(mds)]

cov_rse <- cov_rse[, cov_rse$BrNum %in% rownames(mds)]

#### reorders according to rse_gene$BrNum  (matched the order of rse_gene)
mds = mds[rse_gene$BrNum,1:5]

colData(rse_gene) = cbind(colData(rse_gene), mds)
colData(cov_rse) = cbind(colData(cov_rse), mds)



###################
##### get qSVs
###################

# modJoint = model.matrix(~PrimaryDx + AgeDeath + Sex + mitoRate + rRNA_rate +
# 				totalAssignedGene + RIN, data = colData(rse_gene))
#
### replaced from run_wgcna_combined

modJoint = model.matrix(~PrimaryDx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 +
	mitoRate + rRNA_rate + totalAssignedGene + RIN,
	data=colData(rse_gene))

colnames(modJoint)

#  [1] "(Intercept)"                      "PrimaryDxControl"
#  [3] "PrimaryDxBipolar"                 "BrainRegionsACC"
#  [5] "AgeDeath"                         "SexM"
#  [7] "snpPC1"                           "snpPC2"
#  [9] "snpPC3"                           "mitoRate"
# [11] "rRNA_rate"                        "totalAssignedGene"
# [13] "RIN"                              "PrimaryDxControl:BrainRegionsACC"
# [15] "PrimaryDxBipolar:BrainRegionsACC"

degExprs = log2(assays(cov_rse)$count+1)
k = num.sv(degExprs, modJoint) # 23
print(k) #24
qSV_mat = prcomp(t(degExprs))$x[,1:k]
varExplQsva = jaffelab::getPcaVars(prcomp(t(degExprs)))
varExplQsva[1:k]

#  [1] 67.600  4.820  2.990  2.000  1.420  1.230  1.130  0.805  0.775  0.671
# [11]  0.587  0.528  0.482  0.435  0.351  0.336  0.327  0.296  0.249  0.235
# [21]  0.231  0.225  0.206  0.196


sum(varExplQsva[1:k]) # 87.888%

# [1] 88.125

# model w/o interaction to subset by region
modSep = model.matrix(~PrimaryDx + AgeDeath + Sex  + snpPC1 + snpPC2 + snpPC3 + mitoRate + rRNA_rate +
				totalAssignedGene + RIN, data=colData(rse_gene))

colnames(modSep)

# [1] "(Intercept)"       "PrimaryDxControl"  "PrimaryDxBipolar"
#  [4] "AgeDeath"          "SexM"              "snpPC1"
#  [7] "snpPC2"            "snpPC3"            "mitoRate"
# [10] "rRNA_rate"         "totalAssignedGene" "RIN"

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
### borrows information for
vGene_sACC = voom(dge_sACC,mod_sACC, plot=FALSE)

### limma commads
fitGene_sACC = lmFit(vGene_sACC)
eBGene_sACC = eBayes(fitGene_sACC)
#### need to know  what went into model (modssep) -- case-control and  BP control
#### because of more than 1 component, computing F statistics instead of t=-statistics
outGene_sACC = topTable(eBGene_sACC,coef=2:3,
	p.value = 1,number=nrow(rse_gene))
#### reorderoing based on genes in rse_gene
outGene_sACC = outGene_sACC[rownames(rse_gene),]

## significance levels EXTRACT INDIVIDUAL COMPARISON P-VALUES THAT ARE NOT IN TOP TABLE
pvalMat = as.matrix(eBGene_sACC$p.value)[,2:3]

### check top p-values

head(pvalMat[order(pvalMat[,"PrimaryDxControl"]), ])

qvalMat = pvalMat
qvalMat[,1:2] = p.adjust(pvalMat[,1:2],method="fdr")
colnames(pvalMat) = paste0("P_",colnames(pvalMat))
colnames(qvalMat) = paste0("q_",colnames(qvalMat))

outGene_sACC = cbind(outGene_sACC,cbind(pvalMat, qvalMat))
head(outGene_sACC)
sum(outGene_sACC$q_PrimaryDxControl < 0.05)


# 656
print("Gene_sACC")
sum(outGene_sACC$q_PrimaryDxControl < 0.05)
# [1] 656
sum(outGene_sACC$q_PrimaryDxControl < 0.01)
# [1] 205

sum(outGene_sACC$q_PrimaryDxBipolar < 0.05)
# [1] 218
sum(outGene_sACC$q_PrimaryDxBipolar < 0.01)
# [1] 55


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
# 94
print("Gene_Amyg")
sum(outGene_Amyg$q_PrimaryDxControl < 0.05)
# [1] 94
sum(outGene_Amyg$q_PrimaryDxControl < 0.01)
# [1] 31

sum(outGene_Amyg$q_PrimaryDxBipolar < 0.05)
# [1] 103
sum(outGene_Amyg$q_PrimaryDxBipolar < 0.01)
# [1] 53

save(outGene_Amyg, outGene_sACC, file="qSVA_MDD_gene_DEresults.rda")



#################
### Exon ########
#################

##### sACC ######
dge_sACC = DGEList(counts = assays(rse_exon[,sACC_Index])$counts,
	genes = rowData(rse_exon))
dge_sACC = calcNormFactors(dge_sACC)
vGene_sACC = voom(dge_sACC,mod_sACC, plot=FALSE)

fitGene_sACC = lmFit(vGene_sACC)
eBGene_sACC = eBayes(fitGene_sACC)
outExon_sACC = topTable(eBGene_sACC,coef=2:3,
	p.value = 1,number=nrow(rse_exon))
outExon_sACC = outExon_sACC[rownames(rse_exon),]

## significance levels
pvalMat = as.matrix(eBGene_sACC$p.value)[,2:3]
qvalMat = pvalMat
qvalMat[,1:2] = p.adjust(pvalMat[,1:2],method="fdr")
colnames(pvalMat) = paste0("P_",colnames(pvalMat))
colnames(qvalMat) = paste0("q_",colnames(qvalMat))

outExon_sACC = cbind(outExon_sACC,cbind(pvalMat, qvalMat))
sum(outExon_sACC$q_PrimaryDxControl < 0.05)
# 3776

print("Exon_sACC")
sum(outExon_sACC$q_PrimaryDxControl < 0.05)
# [1] 3776
sum(outExon_sACC$q_PrimaryDxControl < 0.01)
# [1] 938
length(unique(outExon_sACC[which(outExon_sACC$q_PrimaryDxControl < 0.05),]$ensemblID))
# [1] 1358

sum(outExon_sACC$q_PrimaryDxBipolar < 0.05)
# [1] 2045
sum(outExon_sACC$q_PrimaryDxBipolar < 0.01)
# [1] 427
length(unique(outExon_sACC[which(outExon_sACC$q_PrimaryDxBipolar < 0.05),]$ensemblID))
# [1] 844


##### Amygdala ######
dge_Amyg = DGEList(counts = assays(rse_exon[,Amyg_Index])$counts,
	genes = rowData(rse_exon))
dge_Amyg = calcNormFactors(dge_Amyg)
vGene_Amyg = voom(dge_Amyg,mod_Amyg, plot=FALSE)

fitGene_Amyg = lmFit(vGene_Amyg)
eBGene_Amyg = eBayes(fitGene_Amyg)
outExon_Amyg = topTable(eBGene_Amyg,coef=2:3,
	p.value = 1,number=nrow(rse_exon))
outExon_Amyg = outExon_Amyg[rownames(rse_exon),]

## significance levels
pvalMat = as.matrix(eBGene_Amyg$p.value)[,2:3]
qvalMat = pvalMat
qvalMat[,1:2] = p.adjust(pvalMat[,1:2],method="fdr")
colnames(pvalMat) = paste0("P_",colnames(pvalMat))
colnames(qvalMat) = paste0("q_",colnames(qvalMat))

outExon_Amyg = cbind(outExon_Amyg,cbind(pvalMat, qvalMat))
sum(outExon_Amyg$q_PrimaryDxControl < 0.05)
# 1219
print("Exon_Amyg")

sum(outExon_Amyg$q_PrimaryDxControl < 0.05)
# [1] 1219
sum(outExon_Amyg$q_PrimaryDxControl < 0.01)
# [1] 315
length(unique(outExon_Amyg[which(outExon_Amyg$q_PrimaryDxControl < 0.05),]$ensemblID))
# [1] 520

sum(outExon_Amyg$q_PrimaryDxBipolar < 0.05)
# [1] 2102
sum(outExon_Amyg$q_PrimaryDxBipolar < 0.01)
# [1] 612
length(unique(outExon_Amyg[which(outExon_Amyg$q_PrimaryDxBipolar < 0.05),]$ensemblID))
# [1] 837

save(outExon_Amyg, outExon_sACC, file="qSVA_MDD_exon_DEresults.rda")




#################
### Jxn ########
#################

##### sACC ######
dge_sACC = DGEList(counts = assays(rse_jxn[,sACC_Index])$counts,
	genes = rowData(rse_jxn))
dge_sACC = calcNormFactors(dge_sACC)
vGene_sACC = voom(dge_sACC,mod_sACC, plot=FALSE)

fitGene_sACC = lmFit(vGene_sACC)
eBGene_sACC = eBayes(fitGene_sACC)
outJxn_sACC = topTable(eBGene_sACC,coef=2:3,
	p.value = 1,number=nrow(rse_jxn))
outJxn_sACC = outJxn_sACC[rownames(rse_jxn),]

## significance levels
pvalMat = as.matrix(eBGene_sACC$p.value)[,2:3]
qvalMat = pvalMat
qvalMat[,1:2] = p.adjust(pvalMat[,1:2],method="fdr")
colnames(pvalMat) = paste0("P_",colnames(pvalMat))
colnames(qvalMat) = paste0("q_",colnames(qvalMat))

outJxn_sACC = cbind(outJxn_sACC,cbind(pvalMat, qvalMat))
sum(outJxn_sACC$q_PrimaryDxControl < 0.05)
# 1331


##### Amygdala ######
dge_Amyg = DGEList(counts = assays(rse_jxn[,Amyg_Index])$counts,
	genes = rowData(rse_jxn))
dge_Amyg = calcNormFactors(dge_Amyg)
vGene_Amyg = voom(dge_Amyg,mod_Amyg, plot=FALSE)

fitGene_Amyg = lmFit(vGene_Amyg)
eBGene_Amyg = eBayes(fitGene_Amyg)
outJxn_Amyg = topTable(eBGene_Amyg,coef=2:3,
	p.value = 1,number=nrow(rse_jxn))
outJxn_Amyg = outJxn_Amyg[rownames(rse_jxn),]

## significance levels
pvalMat = as.matrix(eBGene_Amyg$p.value)[,2:3]
qvalMat = pvalMat
qvalMat[,1:2] = p.adjust(pvalMat[,1:2],method="fdr")
colnames(pvalMat) = paste0("P_",colnames(pvalMat))
colnames(qvalMat) = paste0("q_",colnames(qvalMat))

outJxn_Amyg = cbind(outJxn_Amyg,cbind(pvalMat, qvalMat))
sum(outJxn_Amyg$q_PrimaryDxControl < 0.05)
# 529


## remove novel
outJxn_sACC_anno = outJxn_sACC[which(outJxn_sACC$Class != "Novel"),]
outJxn_Amyg_anno = outJxn_Amyg[which(outJxn_Amyg$Class != "Novel"),]

print("Jxn_sACC")

## sACC
sum(outJxn_sACC_anno$q_PrimaryDxControl < 0.05)
# [1] 1113
sum(outJxn_sACC_anno$q_PrimaryDxControl < 0.01)
# [1] 487
length(unique(outJxn_sACC_anno[which(outJxn_sACC_anno$q_PrimaryDxControl < 0.05),]$newGeneID))
# [1] 714

sum(outJxn_sACC_anno$q_PrimaryDxBipolar < 0.05)
# [1] 686
sum(outJxn_sACC_anno$q_PrimaryDxBipolar < 0.01)
# [1] 249
length(unique(outJxn_sACC_anno[which(outJxn_sACC_anno$q_PrimaryDxBipolar < 0.05),]$newGeneID))
# [1] 481


print("Jxn_Amyg")

## AMYG
sum(outJxn_Amyg_anno$q_PrimaryDxControl < 0.05)
# [1] 361
sum(outJxn_Amyg_anno$q_PrimaryDxControl < 0.01)
# [1] 121
length(unique(outJxn_Amyg_anno[which(outJxn_Amyg_anno$q_PrimaryDxControl < 0.05),]$newGeneID))
# [1] 258

sum(outJxn_Amyg_anno$q_PrimaryDxBipolar < 0.05)
# [1] 647
sum(outJxn_Amyg_anno$q_PrimaryDxBipolar < 0.01)
# [1] 233
length(unique(outJxn_Amyg_anno[which(outJxn_Amyg_anno$q_PrimaryDxBipolar < 0.05),]$newGeneID))
# [1] 468


save(outJxn_Amyg,outJxn_Amyg_anno, outJxn_sACC,outJxn_sACC_anno, file="qSVA_MDD_jxn_DEresults.rda")




#################
### Tx ########
## we do not use voom  which  required count data
#################

txExprs = log2(assays(rse_tx)$tpm + 1)

##### sACC ######
fitGene_sACC = lmFit(txExprs[,sACC_Index], mod_sACC)
eBGene_sACC = eBayes(fitGene_sACC)
outTx_sACC = topTable(eBGene_sACC,coef=2:3,
	p.value = 1,number=nrow(rse_tx))
outTx_sACC = outTx_sACC[rownames(rse_tx),]

## significance levels
pvalMat = as.matrix(eBGene_sACC$p.value)[,2:3]
qvalMat = pvalMat
qvalMat[,1:2] = p.adjust(pvalMat[,1:2],method="fdr")
colnames(pvalMat) = paste0("P_",colnames(pvalMat))
colnames(qvalMat) = paste0("q_",colnames(qvalMat))

outTx_sACC = cbind(rowData(rse_tx)[,c(5:8)],outTx_sACC,cbind(pvalMat, qvalMat))
sum(outTx_sACC$q_PrimaryDxControl < 0.05)
# 1457

print("Tx_sACC")


sum(outTx_sACC$q_PrimaryDxControl < 0.05)
# [1] 1457
sum(outTx_sACC$q_PrimaryDxControl < 0.01)
# [1] 386
length(unique(outTx_sACC[which(outTx_sACC$q_PrimaryDxControl < 0.05),]$gene_id))
# [1] 1384

sum(outTx_sACC$q_PrimaryDxBipolar < 0.05)
# [1] 783
sum(outTx_sACC$q_PrimaryDxBipolar < 0.01)
# [1] 216
length(unique(outTx_sACC[which(outTx_sACC$q_PrimaryDxBipolar < 0.05),]$gene_id))
# [1] 761

##### Amygdala ######
fitGene_Amyg = lmFit(txExprs[,Amyg_Index], mod_Amyg)
eBGene_Amyg = eBayes(fitGene_Amyg)
outTx_Amyg = topTable(eBGene_Amyg,coef=2:3,
	p.value = 1,number=nrow(rse_tx))
outTx_Amyg = outTx_Amyg[rownames(rse_tx),]

## significance levels
pvalMat = as.matrix(eBGene_Amyg$p.value)[,2:3]
qvalMat = pvalMat
qvalMat[,1:2] = p.adjust(pvalMat[,1:2],method="fdr")
colnames(pvalMat) = paste0("P_",colnames(pvalMat))
colnames(qvalMat) = paste0("q_",colnames(qvalMat))

outTx_Amyg = cbind(rowData(rse_tx)[,c(5:8)],outTx_Amyg,cbind(pvalMat, qvalMat))
sum(outTx_Amyg$q_PrimaryDxControl < 0.05)
# 2101

print("Tx_Amyg")


sum(outTx_Amyg$q_PrimaryDxControl < 0.05)
# [1] 2101
sum(outTx_Amyg$q_PrimaryDxControl < 0.01)
# [1] 521
length(unique(outTx_Amyg[which(outTx_Amyg$q_PrimaryDxControl < 0.05),]$gene_id))
# [1] 1947

sum(outTx_Amyg$q_PrimaryDxBipolar < 0.05)
# [1] 2297
sum(outTx_Amyg$q_PrimaryDxBipolar < 0.01)
# [1] 616
length(unique(outTx_Amyg[which(outTx_Amyg$q_PrimaryDxBipolar < 0.05),]$gene_id))
# [1] 2135

save(outTx_Amyg, outTx_sACC, file="qSVA_MDD_tx_DEresults.rda")



#sgejobs::job_single('qSV_model_DE_analysis', create_shell = TRUE, memory = '80G', command = "Rscript qSV_model_DE_analysis.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
