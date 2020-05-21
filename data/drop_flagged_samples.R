library("sessioninfo")
library("SummarizedExperiment")
library("jaffelab")

## style this script
#styler::style_file("drop_flagged_samples.R", transformers = styler::tidyverse_style(indent_by = 4))

## load rse data
load("rse_gene_GoesZandi_n1093.rda", verbose = TRUE)
load("rse_exon_GoesZandi_n1093.rda", verbose = TRUE)
load("rse_jxn_GoesZandi_n1093.rda", verbose = TRUE)
load("rse_tx_GoesZandi_n1093.rda", verbose = TRUE)

## Find samples to drop
load("../genotype_data/flagged_samples_mddseq.Rdata", verbose = TRUE)

## phenotype data

pd <- colData(rse_gene)
keep_samples <- !pd$RNum %in% checks$SAMPLE_ID
table(keep_samples)
dim(checks)
length(unique(checks$SAMPLE_ID))

### process ERCC for MDD  https://github.com/LieberInstitute/RNAseq-pipeline/blob/master/sh/create_count_objects-human.R#L306-L336


##observed kallisto tpm
sampIDs <- pd$RNum[pd$Experiment=="GoesMDD"]
opt <- list("maindir"="/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data")
erccTPM = sapply(sampIDs, function(x) {
    read.table(file.path(opt$maindir, "Ercc", x, "abundance.tsv"),
               header = TRUE)$tpm
})
rownames(erccTPM) = read.table(file.path(opt$maindir, "Ercc", sampIDs[1], "abundance.tsv"),
                               header = TRUE)$target_id
#check finiteness / change NaNs to 0s
erccTPM[which(is.na(erccTPM), arr.ind = T)] = 0

#expected concentration
spikeIns = read.delim(
    "/users/ajaffe/Lieber/Projects/RNAseq/Ribozero_Compare/ercc_actual_conc.txt",
    as.is = TRUE,
    row.names = 2
)
##match row order
spikeIns = spikeIns[match(rownames(erccTPM), rownames(spikeIns)), ]

pdf(
    file.path('ercc_spikein_check_mix1_drop_flagged_samples.pdf'),
    h = 12,
    w = 18
)
rafalib::mypar(4, 6)
for (i in 1:ncol(erccTPM)) {
    plot(
        log2(10 * spikeIns[, "concentration.in.Mix.1..attomoles.ul."] + 1) ~ log2(erccTPM[, i] +
                                                                                      1),
        xlab = "Kallisto log2(TPM+1)",
        ylab = "Mix 1: log2(10*Concentration+1)",
        main = colnames(erccTPM)[i],
        xlim = c(min(log2(erccTPM + 1)), max(log2(erccTPM + 1)))
    )
    abline(0, 1, lty = 2)
}
dev.off()

mix1conc = matrix(
    rep(spikeIns[, "concentration.in.Mix.1..attomoles.ul."]),
    nc = ncol(erccTPM),
    nr = nrow(erccTPM),
    byrow = FALSE
)
logErr = (log2(erccTPM + 1) - log2(10 * mix1conc + 1))
pd$ERCCsumLogErr <- NA
pd$ERCCsumLogErr[pd$Experiment == "GoesMDD"] = colSums(logErr)

## extract ERCC for BP project samples
metrics <- read.csv("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/read_and_alignment_metrics_zandiHyde_Bipolar_LIBD.csv")
sampIDs_BP <- pd$RNum[pd$Experiment=="ZandiBPD"]
m_bp <- match(sampIDs_BP, metrics$RNum)
table(is.na(m_bp))
pd$ERCCsumLogErr[pd$Experiment == "ZandiBPD"] <- metrics$ERCCsumLogErr[m_bp]

##explore ERCCsumLogErr  (expect no NA)
summary(pd$ERCCsumLogErr)
tapply(pd$ERCCsumLogErr, pd$Experiment, summary)

# $GoesMDD
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -64.9605 -18.6245 -14.9598 -15.8102 -11.8254  -0.4339
#
# $ZandiBPD
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -36.992  -4.586   3.903  10.526  29.208  47.363

## add ERCC
colData(rse_gene) <- pd
colData(rse_exon) <- pd
colData(rse_jxn) <- pd
colData(rse_tx) <- pd



#### subset samples
dim(rse_gene)
rse_gene <- rse_gene[, keep_samples]
dim(rse_gene)

rse_exon <- rse_exon[, keep_samples]
rse_jxn <- rse_jxn[, keep_samples]
rse_tx <- rse_tx[, keep_samples]

## add ancestry
load("../genotype_data/goesHyde_bipolarMdd_Genotypes_n588_mds.rda", verbose = TRUE)
class(mds)
dim(mds)
corner(mds)
table(rse_gene$BrNum %in% rownames(mds))

# FALSE  TRUE
# 12   836
## 6 individual, 12 brain regions absent in mds file

table(rse_gene$BrNum %in% rownames(mds), rse_gene$BrainRegion)
# Amygdala sACC
# FALSE        6    6
# TRUE       415  421

## keep samples with genotypes
rse_gene <- rse_gene[, rse_gene$BrNum %in% rownames(mds)]
cov_rse <- cov_rse[, cov_rse$BrNum %in% rownames(mds)]

dim(rse_gene)
# [1] 25212   836
addmargins(table(rse_gene$Dx, rse_gene$BrainRegion))
# Amygdala sACC Sum
# MDD          234  227 461
# Control      181  194 375
# Sum          415  421 836

#### reorders according to rse_gene$BrNum
mds = mds[rse_gene$BrNum,1:5]
dim(mds)
# [1] 836   5

#colnames(mds) = paste0("snpPC", 1:5)
colData(rse_gene) = cbind(colData(rse_gene), mds)


## save data
dir.create("without_flagged", showWarnings = FALSE)
save(rse_gene, file = file.path("without_flagged", "rse_gene_bp_mdd.Rdata"))
save(rse_exon, file = file.path("without_flagged", "rse_exon_bp_mdd.Rdata"))
save(rse_jxn, file = file.path("without_flagged", "rse_jxn_bp_mdd.Rdata"))
save(rse_tx, file = file.path("without_flagged", "rse_tx_bp_mdd.Rdata"))


# sgejobs::job_single('drop_flagged_samples', create_shell = TRUE, memory = '40G', command = "Rscript drop_flagged_samples.R")


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
