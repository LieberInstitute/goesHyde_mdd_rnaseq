################################## check brain regions
## load packages
library(jaffelab)
library(GenomicRanges)
library(SummarizedExperiment)
library(rtracklayer)
library(RColorBrewer)


## load phenotype and alignment data
pd = read.csv("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/read_and_alignment_metrics_goesHyde_MDD.csv", stringsAsFactors=FALSE, row.names=1)

# One row per BrNum
pdUnique = pd[!duplicated(pd$BrNum),2:6]
pdUnique = pdUnique[order(pdUnique$BrNum),]


# in newest batch?
newBatch = read.csv("190708-190627AS-01_Ran Tao_Omni2.5v1.csv", stringsAsFactors=FALSE)
newBatch = newBatch[!is.na(newBatch$Index),]

pdUnique$inNewBatch = pdUnique$BrNum %in% newBatch$BrNumber


# already in geno?
fam = read.table("/dcl01/lieber/ajaffe/Brain/Imputation/Merged/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_Thrice.fam",
    as.is=TRUE)
colnames(fam) = c("FID", "IID", "MID", "PID", "SEX","PHENO")
fam$BrNum = fam$FID
fam$BrNum[fam$FID == "Omni2pt5"] = fam$IID[fam$FID == "Omni2pt5"]
fam$BrNum = ss(fam$BrNum, "_")

pdUnique$inGeno = pdUnique$BrNum %in% fam$BrNum

write.csv(pdUnique, file="brNums_need_genotype.csv")









