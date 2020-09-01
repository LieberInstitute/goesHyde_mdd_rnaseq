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
qc = read.csv("../qc_checks/qc_dropping_results.csv", stringsAsFactors = FALSE, row.names = 1)
qc = qc[rowSums(qc[,c("dropMetrics","dropRegion","dropRace")])>0,]
# filter from ppt_plot
rse_gene = rse_gene[,-which(rse_gene$RNum %in% rownames(qc))]
#other filters 
rse_gene = rse_gene[,which(!rse_gene$RNum %in% c("R17538","R18853")| rse_gene$PrimaryDx != "Other")] 
message("Drop: ", nrow(qc))

rse_gene$PrimaryDx = factor(rse_gene$PrimaryDx, levels = c("MDD", "Control","Bipolar"))

tempRpkm = recount::getRPKM(rse_gene, "Length")
rowData(rse_gene)$meanExprs = rowMeans(tempRpkm)

pd = colData(rse_gene)

n = ncol(rse_gene)
message("Remaining Samples: ", n)

table(pd$BrainRegion,pd$PrimaryDx)
table(pd$Sex, pd$PrimaryDx)
summary(pd$AgeDeath)

# > table(pd$BrainRegion,pd$PrimaryDx)
# 
# MDD Control Bipolar
# Amygdala 234     188     123
# sACC     230     202     123
# > table(pd$Sex, pd$PrimaryDx)
# 
# MDD Control Bipolar
# F 156      77      96
# M 308     313     150
# > summary(pd$AgeDeath)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 17.37   34.53   47.09   46.55   55.87   95.27 



## Save rse_gene
save(rse_gene, file=paste0("rse_gene_GoesZandi_n",n ,".rda"))


##### Exons

#load objects
load('../preprocessed_data/rse_exon_goesHyde_MDD_n634.Rdata', verbose=TRUE)
rse_mdd = rse_exon	
rse_mdd$Experiment = "psychENCODE_MDD"
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rse_exon_zandiHyde_Bipolar_LIBD_n540.Rdata", verbose=TRUE)
rse_bip = rse_exon		
rse_bip$Experiment = "psychENCODE_BP"

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
m = findMatches(rowRanges(rse_mdd), rowRanges(rse_bip))
rse_mdd = rse_mdd[queryHits(m),]
rse_bip = rse_bip[subjectHits(m),]
rowData(rse_bip) = rowData(rse_mdd)

### combine
rse_both = cbind(rse_mdd, rse_bip)
rowData(rse_both)$meanExprs = NULL
colnames(rse_both) <- paste0(rse_both$RNum,"_",rse_both$Experiment)
rse_exon = rse_both

# drop samples
rse_exon = rse_exon[,rownames(colData(rse_gene))]
tempRpkm = recount::getRPKM(rse_exon, "Length")
rowData(rse_exon)$meanExprs = rowMeans(tempRpkm)

save(rse_exon, file=paste0("rse_exon_GoesZandi_n",n,".rda"))



##### Junctions

#load objects
load('../preprocessed_data/rse_jx_goesHyde_MDD_n634.Rdata', verbose=TRUE)
rse_mdd = rse_jx			
rse_mdd$Experiment = "psychENCODE_MDD"
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rse_jx_zandiHyde_Bipolar_LIBD_n540.Rdata", verbose=TRUE)
rse_bip = rse_jx				
rse_bip$Experiment = "psychENCODE_BP"

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
m = findMatches(rowRanges(rse_mdd), rowRanges(rse_bip))
rse_mdd = rse_mdd[queryHits(m),]
rse_bip = rse_bip[subjectHits(m),]
rowData(rse_bip) = rowData(rse_mdd)

### combine
rse_both = cbind(rse_mdd, rse_bip)
colnames(rse_both) <- paste0(rse_both$RNum,"_",rse_both$Experiment)
rse_jxn = rse_both

# drop samples
rse_jxn = rse_jxn[,rownames(colData(rse_gene))]

rowData(rse_jxn)$Length = 100
tempRpkm = recount::getRPKM(rse_jxn, "Length")
rowData(rse_jxn)$meanExprs = rowMeans(tempRpkm)

save(rse_jxn, file=paste0("rse_jxn_GoesZandi_n",n,".rda"))


##### Transcript

#load objects
load('../preprocessed_data/rse_tx_goesHyde_MDD_n634.Rdata', verbose=TRUE)
rse_mdd = rse_tx	
rse_mdd$Experiment = "psychENCODE_MDD"
# load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rse_tx_zandiHyde_Bipolar_LIBD_n540.Rdata", verbose=TRUE)
# rse_bip = rse_tx
# rse_bip$Experiment = "psychENCODE_BP"

rse_bip = rse_gene[,rse_gene$Experiment == "psychENCODE_BP"]
## also load counts
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rawCounts_zandiHyde_tx.rda", verbose=TRUE)
colnames(txNumReads) = ss(colnames(txNumReads), "_")
txNumReads = txNumReads[,rse_bip$RNum]
colnames(txNumReads) = colnames(rse_bip)
assays(rse_bip, withDimnames=FALSE)$counts <- as.matrix(txNumReads)

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

