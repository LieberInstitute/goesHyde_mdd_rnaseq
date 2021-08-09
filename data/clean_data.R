#####################################

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)
library(edgeR)
library(GenomicRanges)
library(rtracklayer)
library(Matrix)
library(reshape2)
library(data.table)
library(here)

## For styling this script
# styler::style_file("clean_data.R", transformers = biocthis::bioc_style())

# load data
load("rse_gene_raw_GoesZandi.rda", verbose = TRUE)
qc <- read.csv(here("qc_checks","qc_dropping_results.csv"), stringsAsFactors = FALSE, row.names = 1)
message("Samples match: ", all(rse_gene$RNum == qc$RNum))

#Correct BrainRegion
rse_gene$BrainRegion <- qc$BrainRegion

## Add ERCCsumLogErr
rse_gene$ERCCsumLogErr <- qc$ERCCsumLogErr

## Extract rse_gene_ol
rse_gene_ol <- rse_gene[, which(rse_gene$overlap)]
message("Overlapping control samples: ", ncol(rse_gene_ol))

## drop
qc <- qc[qc$dropSum > 0, ]
message("Drop ", nrow(qc), " from qc_checks")
# filter from qc_checks results
rse_gene <- rse_gene[, -which(rse_gene$RNum %in% qc$RNum)]
n <- ncol(rse_gene)
# other filters
other_filter <- which(!rse_gene$RNum %in% c("R17538", "R18853") &
                        rse_gene$PrimaryDx != "Other"&
                        !(rse_gene$overlap & rse_gene$Experiment == "psychENCODE_BP"))
message("Drop ", n-length(other_filter), " for other filters")
rse_gene <- rse_gene[, other_filter]

n <- ncol(rse_gene)
message("Remaining Samples: ", n)

rse_gene$PrimaryDx <- factor(rse_gene$PrimaryDx, levels = c("MDD", "Control", "Bipolar"))
rse_gene_ol$PrimaryDx <- factor(rse_gene_ol$PrimaryDx, levels = c("MDD", "Control", "Bipolar"))

tempRpkm <- recount::getRPKM(rse_gene, "Length")
rowData(rse_gene)$meanExprs <- rowMeans(tempRpkm)

tempRpkm <- recount::getRPKM(rse_gene_ol, "Length")
rowData(rse_gene_ol)$meanExprs <- rowMeans(tempRpkm)
## Confirm bam files exist
# message("All bam files exist: ", all(file.exists(rse_gene$bamFile)))

## Save rse_gene and rse_gene_ol
save(rse_gene, file = "rse_gene_GoesZandi.rda")
save(rse_gene_ol, file = "rse_gene_overlap_GoesZandi.rda")

# create objects from rse_gene
pd <- colData(rse_gene)

table(pd$BrainRegion, pd$PrimaryDx)
table(pd$Sex, pd$PrimaryDx)
summary(pd$AgeDeath)

# > table(pd$BrainRegion, pd$PrimaryDx)
# 
# MDD Control Bipolar
# Amygdala 231     187     122
# sACC     228     200     123
# > table(pd$Sex, pd$PrimaryDx)
# 
# MDD Control Bipolar
# F 155      76      96
# M 304     311     149
# > summary(pd$AgeDeath)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 17.37   34.62   47.21   46.58   55.86   95.27


pd_mdd <- pd[pd$Experiment == "psychENCODE_MDD", ]
pd_bip <- pd[pd$Experiment == "psychENCODE_BP", ]
samples_mdd <- rownames(pd_mdd)
samples_bip <- rownames(pd_bip)

message("Samples mmd: ", length(samples_mdd))
message("Samples bip: ", length(samples_bip))

#### Exons ####

# load objects
load(here("preprocessed_data","rse_exon_goesHyde_MDD_n634.Rdata"), verbose = TRUE)
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

save(rse_exon, file = "rse_exon_GoesZandi.rda")
rm(rse_exon)

#### Junctions ####
# load objects
load(here("preprocessed_data","rse_jx_goesHyde_MDD_n634.Rdata"), verbose = TRUE)
rse_mdd <- rse_jx
rse_mdd$Experiment <- "psychENCODE_MDD"
colnames(rse_mdd) <- paste0(rse_mdd$SAMPLE_ID, "_", rse_mdd$Experiment)

dim(rse_mdd)
# [1] 2272461     634
log10(nrow(rse_mdd) * ncol(rse_mdd))
# [1] 9.158586
print(object.size(rse_mdd),units = "auto")
# 12.3 Gb
class(assays(rse_mdd)$counts)
# [1] "matrix" "array" 

load("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/rse_jx_zandiHyde_Bipolar_LIBD_n540.Rdata", verbose = TRUE)
rse_bip <- rse_jx
rse_bip$Experiment <- "psychENCODE_BP"
colnames(rse_bip) <- paste0(rse_bip$RNum, "_", rse_bip$Experiment)

dim(rse_bip)
# [1] 2314770     540
log10(nrow(rse_bip) * ncol(rse_bip))
# [1] 9.096902
print(object.size(rse_bip),units = "auto")
# 4.8 Gb
class(assays(rse_bip)$counts)
# [1] "S4Vectors"

# filter samples
rse_mdd <- rse_mdd[, samples_mdd]
dim(rse_mdd)
# [1] 2272461     634
rse_bip <- rse_bip[, samples_bip]
dim(rse_bip)
# [1] 2314770     540

## Combine rowRanges to find union of jxns
all_jxn <- unique(c(rowRanges(rse_mdd), rowRanges(rse_bip)))
length(all_jxn)
# [1] 2760639

## Create union matrix
jxn_add <- all_jxn[!names(all_jxn) %in% rownames(rse_mdd),]
jxn_mat_mdd <- rbind(as.matrix(assays(rse_mdd)$counts), 
                     matrix(data = 0, 
                            ncol = length(samples_mdd), 
                            nrow = length(jxn_add),
                            dimnames = list(names(jxn_add),samples_mdd)))
dim(jxn_mat_mdd)
# [1] 2760639     588

jxn_add <- all_jxn[!names(all_jxn) %in% rownames(rse_bip),]
jxn_mat_bip <- as.matrix(rbind(assays(rse_bip)$counts, 
                     matrix(data = 0, 
                            ncol = length(samples_bip), 
                            nrow = length(jxn_add),
                            dimnames = list(names(jxn_add),samples_bip))))
dim(jxn_mat_bip)
# [1] 2760639     503

jxn_mat_both <- cbind(jxn_mat_mdd, jxn_mat_bip)
dim(jxn_mat_both)
class(jxn_mat_both)
jxn_mat_both[1:5,499:503]
pryr::object_size(jxn_mat_both)
# 24.3 GB

## define and populate new SE
rse_jxn <- SummarizedExperiment(assays = list(counts = jxn_mat_both), rowRanges = all_jxn, colData = pd)
rm(jxn_mat_mdd, jxn_mat_bip, jxn_mat_both, rse_bip, rse_mdd, rse_jx)

log10(nrow(rse_jxn) * as.numeric(ncol(rse_jxn)))
pryr::object_size(rse_jxn)
# 24.8 GB
rowData(rse_jxn)$Length <- 100
tempRpkm <- recount::getRPKM(rse_jxn, "Length")
rowData(rse_jxn)$meanExprs <- rowMeans(tempRpkm)

save(rse_jxn, file = "rse_jxn_GoesZandi.rda")
#### Create Sparse Matrix of Data ####
# melt_mdd <- melt(assays(rse_mdd)$counts)
# melt_mdd <- melt_mdd[melt_mdd$value != 0,]
# melt_bip <- melt(as.matrix(as.data.frame(assays(rse_bip)$counts)))
# melt_bip <- melt_bip[melt_bip$value != 0,]
# melt_all <- rbindlist(list(melt_mdd, melt_bip))
# # % non-zero data
# 100*nrow(melt_all)/(as.numeric(length(all_jxn))*nrow(pd))
# # [1] 11.75345
# rm(melt_mdd, melt_bip)
# 
# jxn_names_i <- match(melt_all$Var1, names(all_jxn)) 
# sample_names_j <- match(melt_all$Var2, rownames(pd))
# 
# sparse_jxn <- sparseMatrix(i = jxn_names_i, j = sample_names_j, x = melt_all$value)

##define and populate new SE
# rse_jxn <- SummarizedExperiment(assays = list(counts = sparse_jxn), rowRanges = all_jxn, colData = pd)

# rowData(rse_jxn)$Length <- 100
# tempRpkm <- recount::getRPKM(rse_jxn, "Length")
# rowData(rse_jxn)$meanExprs <- rowMeans(tempRpkm)
# pryr::object_size(rse_jxn)
## save
# save(rse_jxn, file = "rse_jxn_sparse_GoesZandi.rda")
## make test data
# rse_jxn_test <- rse_jxn[1:100000,]
# save(rse_jxn_test, file = "rse_jxn_sparse_test_GoesZandi.rda")
## make non-sparse
# pryr::object_size(rse_jxn)
# 4.65 GB
# sparsity(assays(rse_jxn)$counts)

# assays(rse_jxn)$counts <- as.matrix(assays(rse_jxn)$counts)
# save(rse_jxn, file = "rse_jxn_GoesZandi.rda")
# 
# rm(melt_all, jxn_names_i, sample_names_j, sparse_jxn, all_jxn, rse_jxn, rse_jxn_test)
#### Transcript ####

# load objects
load(here("preprocessed_data","rse_tx_goesHyde_MDD_n634.Rdata"), verbose = TRUE)
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

save(rse_tx, file = "rse_tx_GoesZandi.rda")


# sgejobs::job_single('clean_data', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript clean_data.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

