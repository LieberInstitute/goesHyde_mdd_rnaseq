#####################################

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)
library(edgeR)


# load data
load("rse_gene_raw_GoesZandi_n1140.rda", verbose = TRUE)

## drop
qc = read.csv("../qc_checks/qc_dropping_results.csv", stringsAsFactors = FALSE)
rse_gene = rse_gene[,which(rse_gene$RNum %in% qc$SAMPLE_ID|
                      rse_gene$Experiment == "psychENCODE_BP"| 
                        rse_gene$RNum %in% c("R17538","R18853") |
                        rse_gene$PrimaryDx == "Other" |
                        rse_gene$overallMapRate <0.5) ]

rse_gene$PrimaryDx = factor(rse_gene$PrimaryDx, levels = c("MDD", "Control","Bipolar"))

tempRpkm = recount::getRPKM(rse_gene, "Length")
rowData(rse_gene)$meanExprs = rowMeans(tempRpkm)

pd = colData(rse_gene)

table(pd$BrainRegion,pd$PrimaryDx)
table(pd$Sex, pd$PrimaryDx)
summary(pd$AgeDeath)

# > table(pd$BrainRegion,pd$PrimaryDx)
# 
# MDD Control Bipolar
# Amygdala 237     190     125
# sACC     231     206     128
# > table(pd$BrainRegion,pd$PrimaryDx)
# 
# MDD Control Bipolar
# Amygdala 237     190     125
# sACC     231     206     128
# > table(pd$Sex, pd$PrimaryDx)
# 
# MDD Control Bipolar
# F 158      79      99
# M 310     317     154
# > summary(pd$AgeDeath)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 17.37   34.56   47.09   46.55   55.88   95.27 

n = ncol(rse_gene)
message("Remaining Samples: n", n)

## Save rse_gene
save(rse_gene, file=paste0("rse_gene_GoesZandi_n",n ,".rda"))


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

