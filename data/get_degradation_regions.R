###

library(jaffelab)
library(rtracklayer)
library(recount)
library(recount.bwtool)
library(BiocParallel)
library(SummarizedExperiment)

## For styling this script
# styler::style_file("get_degradation_regions.R", transformers = biocthis::bioc_style())

## use regions from Zandi bip-seq
bed_overall <- import("/dcl01/lieber/ajaffe/lab/degradation_experiments/Joint/bipseq_sACC_Amygdala_RiboZero/bed/sACC_Plus_Amygdala_RibZero_degradation_top1000.bed")
load("rse_gene_GoesZandi.rda", verbose = TRUE)

mdd_samples <- rse_gene$SAMPLE_ID[rse_gene$Experiment == "psychENCODE_MDD"]
bip_samples <- rse_gene$SAMPLE_ID[rse_gene$Experiment == "psychENCODE_BP"]

## designate bigwigs
forwardBw <- c(
    paste0(
        "../preprocessed_data/Coverage/",
        mdd_samples, ".Forward.bw"
    ),
    paste0(
        "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/Coverage/",
        bip_samples, ".Forward.bw"
    )
)
reverseBw <- c(
    paste0(
        "../preprocessed_data/Coverage/",
        mdd_samples, ".Reverse.bw"
    ),
    paste0(
        "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/Coverage/",
        bip_samples, ".Reverse.bw"
    )
)
all(file.exists(c(forwardBw, reverseBw))) # TRUE

names(forwardBw) <- names(reverseBw) <- c(mdd_samples, bip_samples)

## try coverage tool
covForward <- coverage_bwtool(forwardBw, bed_overall,
    strand = "+",
    sumsdir = "degradation_joint", bpparam = MulticoreParam(8)
)
covForward$bigwig_path <- NULL
covForward$bigwig_file <- NULL

covReverse <- coverage_bwtool(reverseBw, bed_overall,
    strand = "-",
    sumsdir = "degradation_joint", bpparam = MulticoreParam(8)
)
covReverse$bigwig_path <- NULL
covReverse$bigwig_file <- NULL

## combine
cov_rse <- rbind(covForward, covReverse)
rownames(cov_rse) <- rowData(cov_rse)$name
cov_rse <- cov_rse[bed_overall$name, ]

## divide by read length
assays(cov_rse)$counts <- assays(cov_rse)$counts / 100 # divide by read length

## make positive
assays(cov_rse)$counts <- abs(assays(cov_rse)$counts)

## Check order matches
stopifnot(all(colnames(cov_rse) == rse_gene$SAMPLE_ID))

## Add coldata
colData(cov_rse) <- colData(rse_gene)

## use RNum_Experiment as as colnames
colnames(cov_rse) <- paste0(cov_rse$RNum, "_", cov_rse$Experiment)

## save
save(cov_rse, file = "degradation_rse_MDDseq_BiPSeq_BothRegions.Rdata")

# sgejobs::job_single('get_degradation_regions', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "get_degradation_regions.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
