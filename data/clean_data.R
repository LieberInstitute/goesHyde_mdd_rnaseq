#####################################

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)
library(edgeR)
library(GenomicRanges)
library(rtracklayer)

## For styling this script
# styler::style_file("clean_data.R", transformers = biocthis::bioc_style())

# load data
load("rse_gene_raw_GoesZandi_n1140.rda", verbose = TRUE)
qc <- read.csv("../qc_checks/qc_dropping_results.csv", stringsAsFactors = FALSE, row.names = 1)
all(rse_gene$RNum == qc$RNum)

#Correct BrainRegion
qc[qc$BrainRegion != rse_gene$BrainRegion,]
#             Experiment  BrNum AgeDeath Sex Race PrimaryDx BrainRegion RIN Plate dropMetrics dropRegion dropRace
# R14179 psychENCODE_MDD Br1469    28.57   M CAUC   Control    Amygdala 7.6     2       FALSE      FALSE    FALSE
# R17496 psychENCODE_MDD Br1469    28.57   M CAUC   Control        sACC 8.0     2       FALSE      FALSE    FALSE
rse_gene$BrainRegion <- qc$BrainRegion

## Add ERCCsumLogErr
rse_gene$ERCCsumLogErr <- qc$ERCCsumLogErr

## drop
qc <- qc[rowSums(qc[, c("dropMetrics", "dropRegion", "dropRace")]) > 0, ]
# filter from ppt_plot
rse_gene <- rse_gene[, -which(rse_gene$RNum %in% rownames(qc))]
# other filters
rse_gene <- rse_gene[, which(!rse_gene$RNum %in% c("R17538", "R18853") | rse_gene$PrimaryDx != "Other")]
message("Drop: ", nrow(qc))

rse_gene$PrimaryDx <- factor(rse_gene$PrimaryDx, levels = c("MDD", "Control", "Bipolar"))

tempRpkm <- recount::getRPKM(rse_gene, "Length")
rowData(rse_gene)$meanExprs <- rowMeans(tempRpkm)

n <- ncol(rse_gene)
message("Remaining Samples: ", n)

## Add bam files to colData
exp_dir <- rep("goesHyde_mdd_rnaseq",n)
exp_dir[rse_gene$Experiment == "psychENCODE_BP"] <- "zandiHyde_bipolar_rnaseq"

bam <- paste0("/dcl01/lieber/ajaffe/lab/",exp_dir,"/preprocessed_data/HISAT2_out/", rse_gene$SAMPLE_ID, "_accepted_hits.sorted.bam")
message("All bam files exist: ", all(file.exists(bam)))
rse_gene$bam_file <- bam

## Save rse_gene
save(rse_gene, file = paste0("rse_gene_GoesZandi_n", n, ".rda"))

# create objects from rse_gene
pd <- colData(rse_gene)

table(pd$BrainRegion, pd$PrimaryDx)
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


pd_mdd <- pd[pd$Experiment == "psychENCODE_MDD", ]
pd_bip <- pd[pd$Experiment == "psychENCODE_BP", ]
samples_mdd <- rownames(pd_mdd)
samples_bip <- rownames(pd_bip)

message("Samples mmd: ", length(samples_mdd))
message("Samples bip: ", length(samples_bip))

##### Exons

# load objects
load("../preprocessed_data/rse_exon_goesHyde_MDD_n634.Rdata", verbose = TRUE)
rse_mdd <- rse_exon
rse_mdd$Experiment <- "psychENCODE_MDD"
colnames(rse_mdd) <- paste0(rse_mdd$SAMPLE_ID, "_", rse_mdd$Experiment)

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rse_exon_zandiHyde_Bipolar_LIBD_n540.Rdata", verbose = TRUE)
rse_bip <- rse_exon
rse_bip$Experiment <- "psychENCODE_BP"
colnames(rse_bip) <- paste0(rse_bip$RNum, "_", rse_bip$Experiment)

# filter samples
rse_mdd <- rse_mdd[, samples_mdd]
rse_bip <- rse_bip[, samples_bip]

# assign pd
colData(rse_mdd) <- pd_mdd
colData(rse_bip) <- pd_bip

# make rowData consistent
m <- findMatches(rowRanges(rse_mdd), rowRanges(rse_bip))
rse_mdd <- rse_mdd[queryHits(m), ]
rse_bip <- rse_bip[subjectHits(m), ]
rowData(rse_bip) <- rowData(rse_mdd)

### combine
rse_exon <- cbind(rse_mdd, rse_bip)
rowData(rse_exon)$meanExprs <- NULL

tempRpkm <- recount::getRPKM(rse_exon, "Length")
rowData(rse_exon)$meanExprs <- rowMeans(tempRpkm)

save(rse_exon, file = paste0("rse_exon_GoesZandi_n", n, ".rda"))
rm(rse_exon)
##### Junctions

# load objects
load("../preprocessed_data/rse_jx_goesHyde_MDD_n634.Rdata", verbose = TRUE)
rse_mdd <- rse_jx
rse_mdd$Experiment <- "psychENCODE_MDD"
colnames(rse_mdd) <- paste0(rse_mdd$SAMPLE_ID, "_", rse_mdd$Experiment)

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rse_jx_zandiHyde_Bipolar_LIBD_n540.Rdata", verbose = TRUE)
rse_bip <- rse_jx
rse_bip$Experiment <- "psychENCODE_BP"
colnames(rse_bip) <- paste0(rse_bip$RNum, "_", rse_bip$Experiment)

# filter samples
rse_mdd <- rse_mdd[, samples_mdd]
rse_bip <- rse_bip[, samples_bip]

# asssign pd
colData(rse_mdd) <- pd_mdd
colData(rse_bip) <- pd_bip

# make rowData consistent
m <- findMatches(rowRanges(rse_mdd), rowRanges(rse_bip))
rse_mdd <- rse_mdd[queryHits(m), ]
rse_bip <- rse_bip[subjectHits(m), ]
rowData(rse_bip) <- rowData(rse_mdd)

## fix class
assays(rse_mdd)$counts = as.matrix(as.data.frame(assays(rse_mdd)$counts ))
## bipseq is a weird SimpleDFrameList
assays(rse_bip) = list(counts = as.matrix(as.data.frame(assays(rse_bip)$counts )))

### combine
rse_jxn <- cbind(rse_mdd, rse_bip)

rowData(rse_jxn)$Length <- 100
tempRpkm <- recount::getRPKM(rse_jxn, "Length")
rowData(rse_jxn)$meanExprs <- rowMeans(tempRpkm)

save(rse_jxn, file = paste0("rse_jxn_GoesZandi_n", n, ".rda"))


##### Transcript

# load objects
load("../preprocessed_data/rse_tx_goesHyde_MDD_n634.Rdata", verbose = TRUE)
rse_mdd <- rse_tx
rse_mdd$Experiment <- "psychENCODE_MDD"
colnames(rse_mdd) <- paste0(rse_mdd$SAMPLE_ID, "_", rse_mdd$Experiment)

rse_bip <- rse_gene[, rse_gene$Experiment == "psychENCODE_BP"]

## transcript
## build rse_tx for bip data (from  zandiHyde_bipolar_rnaseq/data/select_samples.R)
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rpkmCounts_zandiHyde_Bipolar_LIBD_n540.rda", verbose = TRUE)
rm(list = ls()[!ls() %in% c("n", "rse_mdd", "rse_bip", "txTpm", "samples_mdd", "samples_bip", "pd_mdd", "pd_bip")])

gtf <- import("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf")
tx <- gtf[which(gtf$type == "transcript")]
names(tx) <- tx$transcript_id

identical(names(tx), rownames(txTpm))

colnames(txTpm) <- ss(colnames(txTpm), "_")
txTpm <- txTpm[, rse_bip$RNum] ## filter
colnames(txTpm) <- colnames(rse_bip)
rse_bip <- SummarizedExperiment(
    assays = list("tpm" = txTpm),
    colData = colData(rse_bip), rowRanges = tx
)

## also load counts
load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rawCounts_zandiHyde_tx.rda", verbose = T)
colnames(txNumReads) <- ss(colnames(txNumReads), "_")
txNumReads <- txNumReads[, rse_bip$RNum]
colnames(txNumReads) <- colnames(rse_bip)
assays(rse_bip)$counts <- as.matrix(txNumReads)
colnames(rse_bip) <- paste0(rse_bip$RNum, "_", rse_bip$Experiment)

# filter samples
rse_mdd <- rse_mdd[, samples_mdd]
rse_bip <- rse_bip[, samples_bip]

# asssign pd
colData(rse_mdd) <- pd_mdd
colData(rse_bip) <- pd_bip

# make rowData consistent
stopifnot(identical(ranges(rse_bip), ranges(rse_mdd)))

### combine
rse_tx <- cbind(rse_mdd, rse_bip)

save(rse_tx, file = paste0("rse_tx_GoesZandi_n", n, ".rda"))


# sgejobs::job_single('clean_data', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript clean_data.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

