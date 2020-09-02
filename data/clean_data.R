#####################################

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)
library(edgeR)
library(GenomicRanges)
library(rtracklayer)

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

pd_mdd <- pd[pd$Experiment == "psychENCODE_MDD",]
pd_bip <- pd[pd$Experiment == "psychENCODE_BP",]
samples_mdd <- rownames(pd_mdd)
samples_bip <- rownames(pd_bip)

message("Samples mmd: ", length(RNum_mmd))
message("Samples bip: ", length(RNum_bip))

##### Exons

#load objects
load('../preprocessed_data/rse_exon_goesHyde_MDD_n634.Rdata', verbose=TRUE)
rse_mdd = rse_exon	
rse_mdd$Experiment = "psychENCODE_MDD"
colnames(rse_mdd) <- paste0(rse_mdd$SAMPLE_ID,"_",rse_mdd$Experiment)

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rse_exon_zandiHyde_Bipolar_LIBD_n540.Rdata", verbose=TRUE)
rse_bip = rse_exon		
rse_bip$Experiment = "psychENCODE_BP"
colnames(rse_bip) <- paste0(rse_bip$RNum,"_",rse_bip$Experiment)

# filter samples
rse_mdd <- rse_mdd[,samples_mdd]
rse_bip <- rse_bip[,samples_bip]

# asssign pd
colData(rse_mdd) <- pd_mdd
colData(rse_bip) <- pd_bip

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
tempRpkm = recount::getRPKM(rse_exon, "Length")
rowData(rse_exon)$meanExprs = rowMeans(tempRpkm)

save(rse_exon, file=paste0("rse_exon_GoesZandi_n",n,".rda"))



##### Junctions

#load objects
load('../preprocessed_data/rse_jx_goesHyde_MDD_n634.Rdata', verbose=TRUE)
rse_mdd = rse_jx			
rse_mdd$Experiment = "psychENCODE_MDD"
colnames(rse_mdd) <- paste0(rse_mdd$SAMPLE_ID,"_",rse_mdd$Experiment)

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rse_jx_zandiHyde_Bipolar_LIBD_n540.Rdata", verbose=TRUE)
rse_bip = rse_jx				
rse_bip$Experiment = "psychENCODE_BP"
colnames(rse_bip) <- paste0(rse_bip$RNum,"_",rse_bip$Experiment)

# filter samples
rse_mdd <- rse_mdd[,samples_mdd]
rse_bip <- rse_bip[,samples_bip]

# asssign pd
colData(rse_mdd) <- pd_mdd
colData(rse_bip) <- pd_bip

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

save(rse_jxn, file=paste0("rse_jxn_GoesZandi_n",n,".rda"))


##### Transcript

#load objects
load('../preprocessed_data/rse_tx_goesHyde_MDD_n634.Rdata', verbose=TRUE)
rse_mdd = rse_tx	
rse_mdd$Experiment = "psychENCODE_MDD"
colnames(rse_mdd) <- paste0(rse_mdd$SAMPLE_ID,"_",rse_mdd$Experiment)

rse_bip = rse_gene[,rse_gene$Experiment == "psychENCODE_BP"]

## transcript
## build rse_tx for bip data (from  zandiHyde_bipolar_rnaseq/data/select_samples.R)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rpkmCounts_zandiHyde_Bipolar_LIBD_n540.rda", verbose=TRUE)
rm(list = ls()[!ls() %in% c("n","rse_mdd","rse_bip","rse_gene","txTpm")])

gtf = import("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf")
tx = gtf[which(gtf$type == "transcript")]
names(tx) = tx$transcript_id

identical(names(tx), rownames(txTpm))

colnames(txTpm) = ss(colnames(txTpm), "_")
txTpm = txTpm[,rse_bip$RNum]    ## filter
colnames(txTpm) = colnames(rse_bip)
rse_bip = SummarizedExperiment(
  assays = list('tpm' = txTpm),
  colData = colData(rse_bip), rowRanges = tx)

## also load counts
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rawCounts_zandiHyde_tx.rda", verbose=T)
colnames(txNumReads) = ss(colnames(txNumReads), "_")
txNumReads = txNumReads[,rse_bip$RNum]
colnames(txNumReads) = colnames(rse_bip)
assays(rse_bip)$counts = as.matrix(txNumReads)
colnames(rse_bip) <- paste0(rse_bip$RNum,"_",rse_bip$Experiment)

# filter samples
rse_mdd <- rse_mdd[,samples_mdd]
rse_bip <- rse_bip[,samples_bip]

# asssign pd
colData(rse_mdd) <- pd_mdd
colData(rse_bip) <- pd_bip

# make rowData consistent
stopifnot(identical(ranges(rse_bip), ranges(rse_mdd)))

### combine
rse_both = cbind(rse_mdd, rse_bip)
rse_tx = rse_both

# drop samples
rse_tx = rse_tx[,rownames(colData(rse_gene))]

save(rse_tx, file=paste0("rse_tx_GoesZandi_n",n,".rda"))


# sgejobs::job_single('clean_data', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript clean_data.R")

