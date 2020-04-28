library("sessioninfo")
library("SummarizedExperiment")

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

### process ERCC https://github.com/LieberInstitute/RNAseq-pipeline/blob/master/sh/create_count_objects-human.R#L306-L336


##observed kallisto tpm
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
    file.path(opt$maindir, 'Ercc', 'ercc_spikein_check_mix1.pdf'),
    h = 12,
    w = 18
)
mypar(4, 6)
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
metrics$ERCCsumLogErr = colSums(logErr)





#### subset samples
dim(rse_gene)
rse_gene <- rse_gene[, keep_samples]
dim(rse_gene)

rse_exon <- rse_exon[, keep_samples]
rse_jxn <- rse_jxn[, keep_samples]
rse_tx <- rse_tx[, keep_samples]

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
