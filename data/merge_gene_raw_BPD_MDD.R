#####################################

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)
library(edgeR)

setwd('/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/data/')

###############################
##### Load, clean, combine ####
###############################

#load objects
load('../preprocessed_data/rse_gene_goesHyde_MDD_n634.Rdata', verbose=TRUE)
rse_mdd = rse_gene					# 634
rse_mdd$Experiment = "GoesMDD"
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseGene_n511.rda", verbose=TRUE)
rse_bip = rse_gene					# 511
rse_bip$Experiment = "ZandiBPD"
## 9 Control samples also got used in Goes -> drop from Zandi so there's no duplicates
rse_bip = rse_bip[,-which(rse_bip$RNum %in% colnames(rse_mdd))]	# 502

## combine
# make colData consistent
rse_mdd$AgeDeath = rse_mdd$Age
rse_mdd$RNum = rse_mdd$SAMPLE_ID
rse_mdd$BrNum = as.character(rse_mdd$BrNum)
colKeep = c("RNum","BrNum","Sex","Race","AgeDeath","BrainRegion","RIN","PrimaryDx","overallMapRate","totalAssignedGene","mitoRate","rRNA_rate","Experiment")
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
qc = rbind(qc, "R17779")	# Tiny # of reads
rse_both = rse_both[,-which(rse_both$RNum %in% qc$SAMPLE_ID | rse_both$RNum %in% c("R17538","R18853") |
						rse_both$PrimaryDx == "Other" |
						rse_both$overallMapRate <0.5) ]

rse_both$PrimaryDx = droplevels(rse_both$PrimaryDx)
rse_both$PrimaryDx = relevel(rse_both$PrimaryDx, ref="MDD")

tempRpkm = recount::getRPKM(rse_both, "Length")
rowData(rse_both)$meanExprs = rowMeans(tempRpkm)

pd = colData(rse_both)
table(pd$BrainRegion,pd$PrimaryDx)
           # MDD Control Bipolar
  # Amygdala 235     186     120
  # sACC     228     199     125

table(pd$Sex, pd$PrimaryDx)
    # MDD Control Bipolar
  # F 155      78      96
  # M 308     307     149


summary(pd$AgeDeath)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 17.37   34.44   47.21   46.57   55.86   95.27

rse_gene = rse_both
# save(rse_gene, file="rse_gene_GoesZandi_n1093.rda")


##### Exons

#load objects
load('../preprocessed_data/rse_exon_goesHyde_MDD_n634.Rdata', verbose=TRUE)
rse_mdd = rse_exon					# 634
rse_mdd$Experiment = "GoesMDD"
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseExon_n511.rda", verbose=TRUE)
rse_bip = rse_exon					# 511
rse_bip$Experiment = "ZandiBPD"

## combine
# make colData consistent
rse_mdd$AgeDeath = rse_mdd$Age
rse_mdd$RNum = rse_mdd$SAMPLE_ID
rse_mdd$BrNum = as.character(rse_mdd$BrNum)
colKeep = c("RNum","BrNum","Sex","Race","AgeDeath","BrainRegion","RIN","PrimaryDx","overallMapRate","totalAssignedGene","mitoRate","rRNA_rate","Experiment")
colData(rse_mdd) = colData(rse_mdd)[,colKeep]
colData(rse_bip) = colData(rse_bip)[,colKeep]

# make rowData consistent
m = findMatches(rowRanges(rse_mdd), rowRanges(rse_bip))
rse_mdd = rse_mdd[queryHits(m),]
rse_bip = rse_bip[subjectHits(m),]
rowData(rse_bip) = rowData(rse_mdd)

### combine
rse_both = cbind(rse_mdd, rse_bip)
rowData(rse_both)$meanExprs = NULL
rse_exon = rse_both

# drop samples
rse_exon = rse_exon[,rownames(colData(rse_gene))]

tempRpkm = recount::getRPKM(rse_exon, "Length")
rowData(rse_exon)$meanExprs = rowMeans(tempRpkm)

# save(rse_exon, file="rse_exon_GoesZandi_n1093.rda")



##### Junctions

#load objects
load('../preprocessed_data/rse_jx_goesHyde_MDD_n634.Rdata', verbose=TRUE)
rse_mdd = rse_jx					# 634
rse_mdd$Experiment = "GoesMDD"
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseJxn_n511.rda", verbose=TRUE)
rse_bip = rse_jxn					# 511
rse_bip$Experiment = "ZandiBPD"

## combine
# make colData consistent
rse_mdd$AgeDeath = rse_mdd$Age
rse_mdd$RNum = rse_mdd$SAMPLE_ID
rse_mdd$BrNum = as.character(rse_mdd$BrNum)
colKeep = c("RNum","BrNum","Sex","Race","AgeDeath","BrainRegion","RIN","PrimaryDx","overallMapRate","totalAssignedGene","mitoRate","rRNA_rate","Experiment")
colData(rse_mdd) = colData(rse_mdd)[,colKeep]
colData(rse_bip) = colData(rse_bip)[,colKeep]

# make rowData consistent
m = findMatches(rowRanges(rse_mdd), rowRanges(rse_bip))
rse_mdd = rse_mdd[queryHits(m),]
rse_bip = rse_bip[subjectHits(m),]
rowData(rse_bip) = rowData(rse_mdd)

### combine
rse_both = cbind(rse_mdd, rse_bip)
rse_jxn = rse_both

# drop samples
rse_jxn = rse_jxn[,rownames(colData(rse_gene))]

rowData(rse_jxn)$Length = 100
tempRpkm = recount::getRPKM(rse_jxn, "Length")
rowData(rse_jxn)$meanExprs = rowMeans(tempRpkm)

# save(rse_jxn, file="rse_jxn_GoesZandi_n1093.rda")


##### Transcript

#load objects
load('../preprocessed_data/rse_tx_goesHyde_MDD_n634.Rdata', verbose=TRUE)
rse_mdd = rse_tx					# 634
rse_mdd$Experiment = "GoesMDD"
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseTx_n511.rda", verbose=TRUE)
rse_bip = rse_tx					# 511
rse_bip$Experiment = "ZandiBPD"

## also load counts
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rawCounts_zandiHyde_tx.rda", verbose=T)
colnames(txNumReads) = ss(colnames(txNumReads), "_")
txNumReads = txNumReads[,rse_bip$RNum]
colnames(txNumReads) = colnames(rse_bip)
assays(rse_bip)$counts = as.matrix(txNumReads)

## combine
# make colData consistent
rse_mdd$AgeDeath = rse_mdd$Age
rse_mdd$RNum = rse_mdd$SAMPLE_ID
rse_mdd$BrNum = as.character(rse_mdd$BrNum)
colKeep = c("RNum","BrNum","Sex","Race","AgeDeath","BrainRegion","RIN","PrimaryDx","overallMapRate","totalAssignedGene","mitoRate","rRNA_rate","Experiment")
colData(rse_mdd) = colData(rse_mdd)[,colKeep]
colData(rse_bip) = colData(rse_bip)[,colKeep]

# make rowData consistent
stopifnot(identical(ranges(rse_bip), ranges(rse_mdd)))

### combine
rse_both = cbind(rse_mdd, rse_bip)
rse_tx = rse_both

# drop samples
rse_tx = rse_tx[,rownames(colData(rse_gene))]

# save(rse_tx, file="rse_tx_GoesZandi_n1093.rda")


# sgejobs::job_single('clean_data', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript clean_data.R")

