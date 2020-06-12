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
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rse_gene_zandiHyde_Bipolar_LIBD_n540.Rdata", verbose=TRUE) #changed to 540
rse_bip = rse_gene					# 540
rse_bip$Experiment = "ZandiBPD"

## combine
# make colData consistent
rse_mdd$AgeDeath = rse_mdd$Age
rse_mdd$RNum = rse_mdd$SAMPLE_ID
rse_mdd$BrNum = as.character(rse_mdd$BrNum)
rse_bip$BrainRegion = rse_bip$Brain.Region
colKeep = c("RNum","BrNum","Sex","Race","AgeDeath","BrainRegion","RIN","PrimaryDx","overallMapRate","totalAssignedGene","mitoRate","rRNA_rate","Experiment")
colData(rse_mdd) = colData(rse_mdd)[,colKeep]
colData(rse_bip) = colData(rse_bip)[,colKeep]

# make rowData consistent
rowData(rse_bip)$Symbol = rowData(rse_mdd)$Symbol 	# fill in blank symbols
rowData(rse_bip)$meanExprs = rowData(rse_bip)$gencodeTx = NULL
rowData(rse_mdd)$meanExprs = rowData(rse_mdd)$gencodeTx = NULL

### combine
rse_both = cbind(rse_mdd, rse_bip) #1174

## drop
# qc = read.csv("../qc_checks/qc_dropping_results.csv", stringsAsFactors = FALSE)
# qc = qc[rowSums(qc[,13:16])>0,]
#qc = rbind(qc, "R17779")	# Tiny # of reads
#filter to 1130 samples
# rse_both = rse_both[,-which(rse_both$RNum %in% qc$SAMPLE_ID | #38
#                         rse_both$RNum %in% c("R17538","R18853") | #2
# 						rse_both$PrimaryDx == "Other" | #2
# 						rse_both$overallMapRate <0.5) ] #6

rse_both$PrimaryDx = droplevels(rse_both$PrimaryDx)
rse_both$PrimaryDx = relevel(rse_both$PrimaryDx, ref="MDD")

tempRpkm = recount::getRPKM(rse_both, "Length")
rowData(rse_both)$meanExprs = rowMeans(tempRpkm)

pd = colData(rse_both)
table(pd$BrainRegion,pd$PrimaryDx)

rse_gene = rse_both
save(rse_gene, file="rse_gene_raw_GoesZandi_n1174.rda")

#old filtering
           # MDD Control Bipolar
  # Amygdala 235     186     120
  # sACC     228     199     125
# New filtering
  #          MDD Control Bipolar
  # Amygdala 236     197     125
  # sACC     228     215     129
table(pd$Sex, pd$PrimaryDx)
#old
    # MDD Control Bipolar
  # F 155      78      96
  # M 308     307     149
#new
  #   MDD Control Bipolar
  # F 155      81     101
  # M 309     331     153

summary(pd$AgeDeath)
#old
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 17.37   34.44   47.21   46.57   55.86   95.27
#new
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 17.37   34.47   47.30   46.58   55.88   95.27



# sgejobs::job_single('merge_gene_raw_BPD_MDD', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "merge_gene_raw_BPD_MDD.R")

