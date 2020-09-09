###

library(jaffelab)
library(rtracklayer)
library(recount)
library(recount.bwtool)
library(BiocParallel)
library(SummarizedExperiment)

## use regions from Zandi bip-seq
bed_overall = import("/dcl01/lieber/ajaffe/lab/degradation_experiments/Joint/bipseq_sACC_Amygdala_RiboZero/bed/sACC_Plus_Amygdala_RibZero_degradation_top1000.bed")

load("preprocessed_data/rse_gene_goesHyde_MDD_n634.Rdata", verbose=TRUE)

## 2 samples are getting coverage errors
dropInd = which(rse_gene$SAMPLE_ID %in% c("R17538","R18853"))
rse_gene = rse_gene[,-dropInd]

#######################
##### Combined ########
#######################

## designate bigwigs
forwardBw = paste0("preprocessed_data/Coverage/",
	rse_gene$SAMPLE_ID,".Forward.bw")
reverseBw = paste0("preprocessed_data/Coverage/",
	rse_gene$SAMPLE_ID, ".Reverse.bw")
all(file.exists(c(forwardBw,reverseBw))) # TRUE

names(forwardBw) = names(reverseBw) = rse_gene$SAMPLE_ID

## try coverage tool
covForward = coverage_bwtool(forwardBw, bed_overall, strand = "+", 
	sumsdir = "degradation_joint", bpparam = MulticoreParam(8))
covForward$bigwig_path = NULL
covForward$bigwig_file = NULL

covReverse = coverage_bwtool(reverseBw, bed_overall, strand = "-", 
	sumsdir = "degradation_joint", bpparam = MulticoreParam(8))
covReverse$bigwig_path = NULL
covReverse$bigwig_file = NULL

## combine
cov_rse = rbind(covForward, covReverse)	
rownames(cov_rse) = rowData(cov_rse)$name
cov_rse = cov_rse[bed_overall$name,]

## divide by read length
assays(cov_rse)$counts = assays(cov_rse)$counts/100 # divide by read length

## make positive
assays(cov_rse)$counts = abs(assays(cov_rse)$counts) 

## save to final people
colData(cov_rse) = colData(rse_gene)
save(cov_rse, file = "data/degradation_rse_MDDseq_BothRegions.rda")






