##

library(jaffelab)
library(SummarizedExperiment)
library(recount)
library(readxl)
library(RColorBrewer)
library(sessioninfo)
library(dplyr)
library(limma)
library(rtracklayer)
library(genefilter)

## For styling this script
# styler::style_file("qc_checks.R", transformers = biocthis::bioc_style())

## load phenotype and alignment data
load(
    "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/data/rse_gene_raw_GoesZandi.rda",
    verbose = TRUE
)
pd <- colData(rse_gene) %>% as.data.frame()

metricCols <-c("RNum", "Experiment", "numReads", "Plate", "ERCCsumLogErr")
pd_mdd <- read.csv(
    "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data/read_and_alignment_metrics_goesHyde_MDD.csv",
    stringsAsFactors = FALSE,
    row.names = 1) %>%
    mutate(Experiment = "psychENCODE_MDD",
           ERCCsumLogErr = NA) %>%
    rename(RNum = SAMPLE_ID) %>%
    select(all_of(metricCols))

pd_bp <- read.csv(
    "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/preprocessed_data/read_and_alignment_metrics_zandiHyde_Bipolar_LIBD.csv",
    stringsAsFactors = FALSE,
    row.names = 1) %>%
    mutate(Experiment = "psychENCODE_BP",
           Plate = NA) %>%
    select(all_of(metricCols))

pd_both <- rbind(pd_mdd, pd_bp)

pd <- left_join(pd, pd_both, by = c("RNum", "Experiment"))
rownames(pd) <- pd$RNum

#### ERCC data ####
### https://github.com/LieberInstitute/RNAseq-pipeline/blob/master/sh/create_count_objects-human.R#L306-L336

##observed kallisto tpm
sampIDs_mdd <- pd$RNum[pd$Experiment == "psychENCODE_MDD"]
opt <- list("maindir" = "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/preprocessed_data")
erccTPM = sapply(sampIDs_mdd, function(x) {
    read.table(file.path(opt$maindir, "Ercc", x, "abundance.tsv"),
               header = TRUE)$tpm
})
rownames(erccTPM) = read.table(file.path(opt$maindir, "Ercc", sampIDs_mdd[1], "abundance.tsv"),
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
spikeIns = spikeIns[match(rownames(erccTPM), rownames(spikeIns)),]

mix1conc = matrix(
    rep(spikeIns[, "concentration.in.Mix.1..attomoles.ul."]),
    nc = ncol(erccTPM),
    nr = nrow(erccTPM),
    byrow = FALSE
)
logErr = (log2(erccTPM + 1) - log2(10 * mix1conc + 1))
pd$ERCCsumLogErr[pd$Experiment == "psychENCODE_MDD"] = colSums(logErr)

colnames(erccTPM) <- paste0(colnames(erccTPM),": ",round(colSums(logErr),2))
pdf(
    file.path('pdfs/ercc_spikein_check_mix1_drop_flagged_samples.pdf'),
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

pdf("pdfs/ERCCsumLogErr_boxplot.pdf", h = 5, w = 5)
boxplot(
    pd$ERCCsumLogErr ~ pd$Experiment,
    las = 3,
    xlab = "Experiment",
    outline = FALSE,
    ylab = "ERCC sum LogErr"
)
dev.off() 

summary(pd$ERCCsumLogErr)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -64.961 -15.748 -10.183  -3.817   1.050  47.363

tapply(pd$ERCCsumLogErr, pd$Experiment, summary)
# $psychENCODE_BP
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -36.992  -4.592   4.522  10.706  29.660  47.363 
# 
# $psychENCODE_MDD
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -64.9605 -18.8250 -14.8214 -15.8694 -11.8292  -0.4339 


#### Metrics ####
table(pd$BrainRegion, pd$PrimaryDx)
table(pd$Sex, pd$PrimaryDx)
table(pd$Race, pd$PrimaryDx)
table(table(pd$BrNum))
summary(pd$Age)

# > table(pd$BrainRegion, pd$PrimaryDx)

# Bipolar Control MDD
# Amygdala     125     194 240
# sACC         128     214 239
# > table(pd$Sex, pd$PrimaryDx)
#
# Bipolar Control MDD
# F      99      86 159
# M     154     322 320
# > table(pd$Race, pd$PrimaryDx)
#
# Bipolar Control MDD
# AA         0       5   0
# AS         0       2   0
# CAUC     253     397 479
# HISP       0       4   0
# > table(table(pd$BrNum))
#
# 1   2   3
# 71 524   7
# > summary(pd$Age)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 17.37   34.41   47.15   46.55   55.91   95.27


summary(pd$numReads)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 10284026  87697208 103523346 122609239 138632467 453612210


aL <- min(pd$totalAssignedGene)
aH <- max(pd$totalAssignedGene)
mL <- min(pd$overallMapRate)
mH <- max(pd$overallMapRate)
rnaL <- min(pd$rRNA_rate)
rnaH <- max(pd$rRNA_rate)
mitoL <- min(pd$mitoRate)
mitoH <- max(pd$mitoRate)
rL <- min(pd$numReads)
rH <- max(pd$numReads)

# drop samples?
pdf("pdfs/RIN_check_predrop.pdf", h = 10, w = 10)
par(
    mfcol = c(2, 2),
    mar = c(5, 6, 2, 2),
    cex.axis = 1.8,
    cex.lab = 1.8
)
plot(pd$RIN, pd$overallMapRate, pch = 21, bg = "grey")
abline(h = 0.5, lty = 2)
legend("bottomright", "a)", bty = "n", cex = 2)
plot(pd$RIN, pd$totalAssignedGene, pch = 21, bg = "grey")
abline(h = 0.3, lty = 2)
legend("bottomright", "b)", bty = "n", cex = 2)
plot(pd$rRNA_rate,
     pd$overallMapRate,
     pch = 21,
     bg = "grey")
abline(h = 0.5, lty = 2)
abline(v = 5e-4, lty = 2)
legend("bottomright", "c)", bty = "n", cex = 2)
plot(pd$RIN, log10(pd$numReads), pch = 21, bg = "grey")
abline(h = log10(9e6), lty = 2)
legend("bottomright", "d)", bty = "n", cex = 2)
dev.off()


####################################

## drop samples
pd$dropMetrics <- FALSE
pd$dropMetrics[pd$overallMapRate < 0.5 |
                   pd$totalAssignedGene < .3 | pd$numReads < (10^7.25)] <- TRUE

table(pd$dropMetrics)
# FALSE  TRUE
# 1133     7

####################################

pdf("pdfs/RIN_check_postdrop.pdf", h = 10, w = 10)
par(
    mfcol = c(2, 2),
    mar = c(5, 6, 2, 2),
    cex.axis = 1.8,
    cex.lab = 1.8
)
plot(
    pd$RIN,
    pd$overallMapRate,
    pch = 21,
    bg = pd$dropMetrics + 1,
    cex = pd$dropMetrics + 1,
    ylim = c(mL, mH)
)
abline(h = 0.5, lty = 2)
plot(
    pd$RIN,
    pd$totalAssignedGene,
    pch = 21,
    bg = pd$dropMetrics + 1,
    cex = pd$dropMetrics + 1,
    ylim = c(aL, aH)
)
abline(h = 0.3, lty = 2)
plot(
    pd$rRNA_rate,
    pd$overallMapRate,
    pch = 21,
    bg = pd$dropMetrics + 1,
    cex = pd$dropMetrics + 1,
    xlim = c(rnaL, rnaH),
    ylim = c(mL, mH)
)
abline(h = 0.5, lty = 2)
plot(
    pd$RIN,
    log10(pd$numReads),
    pch = 21,
    bg = pd$dropMetrics + 1,
    cex = pd$dropMetrics + 1,
    ylim = c(log10(rL), log10(rH))
)
abline(h = log10(1e7), lty = 2)
dev.off()


#### ERCC vs. Metrics ####
pdf("pdfs/ERCC_check_postdrop.pdf", h = 10, w = 10)
rafalib::mypar(2,2)
# par(
#     mfcol = c(2, 2),
#     mar = c(5, 6, 2, 2),
#     cex.axis = 1.8,
#     cex.lab = 1.8
# )
plot(
    pd$ERCCsumLogErr,
    pd$overallMapRate,
    pch = 21,
    bg = pd$dropMetrics + 1,
    cex = pd$dropMetrics + 1,
    ylim = c(mL, mH)
)
abline(h = 0.5, lty = 2)
plot(
    pd$ERCCsumLogErr,
    pd$totalAssignedGene,
    pch = 21,
    bg = pd$dropMetrics + 1,
    cex = pd$dropMetrics + 1,
    ylim = c(aL, aH)
)
abline(h = 0.3, lty = 2)
plot(
    pd$ERCCsumLogErr,
    log10(pd$numReads),
    pch = 21,
    bg = pd$dropMetrics + 1,
    cex = pd$dropMetrics + 1,
    ylim = c(log10(rL), log10(rH))
)
abline(h = log10(1e7), lty = 2)
plot(
    pd$ERCCsumLogErr,
    pd$mitoRate,
    pch = 21,
    bg = pd$dropMetrics + 1,
    cex = pd$dropMetrics + 1,
    ylim = c(mitoL, mitoH)
)
# abline(h = log10(1e7), lty = 2)
plot(
    pd$ERCCsumLogErr,
    pd$rRNA_rate,
    pch = 21,
    bg = pd$dropMetrics + 1,
    cex = pd$dropMetrics + 1,
    ylim = c(rnaL, rnaH)
)
# abline(h = 0.5, lty = 2)
dev.off()

## metrics by flowcell

table(pd$Plate) # Plate 2 are the 96 rerun samples

### plot
pdf("pdfs/metrics_by_plate.pdf", h = 5, w = 5)
par(
    mar = c(7, 6, 2, 2),
    cex.axis = 1,
    cex.lab = 1.5,
    cex.main = 2
)
palette(brewer.pal(8, "Dark2"))
boxplot(
    pd$overallMapRate ~ pd$Plate,
    las = 3,
    ylim = c(0, 1),
    xlab = "Plate",
    outline = FALSE,
    ylab = "Overall Map Rate"
)
points(
    pd$overallMapRate ~ jitter(as.numeric(factor(pd$Plate)),
                               amount = 0.15),
    pch = 21,
    bg = as.numeric(factor(pd$BrainRegion))
)
legend(
    "bottomleft",
    levels(factor(pd$BrainRegion)),
    pch = 15,
    col = 1:2,
    cex = .8
)
abline(h = .5, lty = 2, col = "grey")

boxplot(
    pd$totalAssignedGene ~ pd$Plate,
    las = 3,
    ylim = c(0, .7),
    xlab = "Plate",
    outline = FALSE,
    ylab = "Gene Assignement Rate"
)
points(
    pd$totalAssignedGene ~ jitter(as.numeric(factor(pd$Plate)),
                                  amount = 0.15),
    pch = 21,
    bg = as.numeric(factor(pd$BrainRegion))
)
abline(h = .3, lty = 2, col = "grey")

dev.off()


#################################
### regional labeling

gRpkm <- recount::getRPKM(rse_gene, "Length")

## filter low expressing genes
gRpkm <- gRpkm[which(rowMeans(gRpkm) > 0.1),]
yExprs <- log2(gRpkm + 1)

### top 1000 genes different between regions
ngenes <- 100

p <- rowttests(yExprs, as.factor(pd$BrainRegion))
pOrd <- p[order(p$p.value)[1:ngenes],]
ind1000 <- which(rownames(p) %in% rownames(pOrd))
yExprs1000 <- yExprs[ind1000,]

### estimate brain region
amyg <- ifelse(pd$BrainRegion == "Amygdala", 1, 0)
sacc <- ifelse(pd$BrainRegion == "sACC", 1, 0)
mod <- data.frame(model.matrix( ~ amyg + sacc - 1))

fit <- lmFit(yExprs1000, mod)
Xmat <- fit$coef
Dmat <- t(Xmat) %*% Xmat
guess <-
    apply(yExprs1000, 2, function(x)
        solve(Dmat, t(Xmat) %*% x))[2,]
stopifnot(identical(names(guess), paste0(rownames(pd), "_", pd$Experiment)))
pd$guess <- guess

### plot
pdf("pdfs/region_check_100.pdf", h = 6, w = 6)
par(
    mar = c(8, 6, 2, 2),
    cex.axis = 1.5,
    cex.lab = 1.5,
    cex.main = 2
)
palette(brewer.pal(8, "Dark2"))
boxplot(
    pd$guess ~ pd$BrainRegion,
    las = 3,
    ylim = range(pd$guess),
    outline = FALSE,
    ylab = "sACC Identity",
    xlab = ""
)
points(
    pd$guess ~ jitter(as.numeric(factor(pd$BrainRegion)),
                      amount = 0.25),
    pch = 21,
    bg = as.numeric(as.factor(pd$BrainRegion))
)
segments(0, 0.9, 1.5, 0.9, lty = 2, col = "grey")
segments(1.5, 0.3, 2.5, 0.3, lty = 2, col = "grey")
dev.off()

info_cols <-
    c(
        "BrNum",
        "Experiment",
        "AgeDeath",
        "Sex",
        "Race",
        "PrimaryDx",
        "BrainRegion",
        "dropMetrics",
        "guess"
    )
pd[which(pd$BrainRegion == "sACC" & pd$guess < .3), info_cols]
pd[which(pd$BrainRegion == "Amygdala" & pd$guess > 0.9), info_cols]


# > pd[which(pd$BrainRegion=="sACC" & pd$guess<.3),info_cols]
# BrNum      Experiment AgeDeath Sex Race PrimaryDx BrainRegion dropMetrics        guess
# R14179 Br1469 psychENCODE_MDD 28.57000   M CAUC   Control        sACC       FALSE -0.021303859
# R17520 Br1635 psychENCODE_MDD 52.31000   M CAUC       MDD        sACC       FALSE  0.231191685
# R17527 Br1675 psychENCODE_MDD 32.10000   M CAUC       MDD        sACC       FALSE -0.056853078
# R17547 Br1754 psychENCODE_MDD 25.98000   F CAUC       MDD        sACC       FALSE  0.018389137
# R17579 Br6099 psychENCODE_MDD 26.66381   M CAUC       MDD        sACC       FALSE  0.153008741
# R17894 Br8017 psychENCODE_MDD 49.29500   M CAUC       MDD        sACC       FALSE  0.261345113
# R17941 Br8133 psychENCODE_MDD 62.36813   M CAUC       MDD        sACC       FALSE -0.104237322
# R17951 Br2333 psychENCODE_MDD 62.99000   F CAUC   Control        sACC       FALSE -0.063507534
# R17952 Br5549 psychENCODE_MDD 67.47000   M CAUC       MDD        sACC       FALSE  0.008696151
# R18430 Br3863 psychENCODE_MDD 69.33904   M CAUC   Control        sACC       FALSE -0.122722764
# R18458 Br8313 psychENCODE_MDD 55.41672   M CAUC   Control        sACC       FALSE  0.134877418
# R19186 Br1475 psychENCODE_MDD 43.30000   M CAUC       MDD        sACC       FALSE  0.066271916
# R14110 Br1562  psychENCODE_BP 36.15000   M CAUC   Bipolar        sACC       FALSE  0.281349529
# R14175 Br1752  psychENCODE_BP 39.20000   F CAUC   Bipolar        sACC       FALSE  0.061872883
# R14182 Br1558  psychENCODE_BP 23.43000   F CAUC   Control        sACC       FALSE  0.266600470
# R14228 Br2333  psychENCODE_BP 62.99000   F CAUC   Control        sACC       FALSE  0.035792461
# R14235 Br1661  psychENCODE_BP 51.15000   M CAUC   Bipolar        sACC       FALSE  0.274285966
# R14269 Br5555  psychENCODE_BP 21.34000   M CAUC   Control        sACC       FALSE  0.134646892
# R14295 Br2071  psychENCODE_BP 40.62000   F CAUC   Bipolar        sACC       FALSE  0.128219700
# > pd[which(pd$BrainRegion=="Amygdala" & pd$guess>0.9),info_cols]
# BrNum      Experiment AgeDeath Sex Race PrimaryDx BrainRegion dropMetrics    guess
# R17496 Br1469 psychENCODE_MDD 28.57000   M CAUC   Control    Amygdala       FALSE 1.075980
# R18938 Br6293 psychENCODE_MDD 27.85763   F CAUC       MDD    Amygdala       FALSE 0.903392
# R14081 Br5611  psychENCODE_BP 63.82000   F CAUC   Bipolar    Amygdala       FALSE 1.010363
# R15072 Br5939  psychENCODE_BP 55.90965   M CAUC   Bipolar    Amygdala       FALSE 1.058777

## Br1469 labels got switched
pd["R14179", "BrainRegion"] <- "Amygdala"
pd["R17496", "BrainRegion"] <- "sACC"

## Drop the rest
dropInd <- which((pd$BrainRegion == "sACC" & pd$guess < .3) |
                     (pd$BrainRegion == "Amygdala" & pd$guess > 0.9))
pd$dropRegion <- FALSE
pd$dropRegion[dropInd] <- TRUE

table(pd$dropRegion)
# FALSE  TRUE
# 1119    21

######################################################## .
gia <- read.csv("../synapse/genotypeInferredAncestry.csv")
pd$geneticRace <- pd$Race
pd$geneticRace[match(gia$BrNum, pd$BrNum)] <-
    gia$genotypeInferredAncestry

pd$dropRace <- pd$geneticRace != "CAUC"

table(pd$Race, pd$geneticRace)
# AA ambiguous   AS CAUC HISP
# AA      5         0    0    0    0
# AS      0         0    2    0    0
# CAUC    0         2    0 1127    0
# HISP    0         0    0    1    3

table(pd$dropRace)
# FALSE  TRUE
# 1128    12


######################################################## .
info_cols <-
    c(  "RNum",
        "BrNum",
        "Experiment",
        "AgeDeath",
        "Sex",
        "Race",
        "PrimaryDx",
        "BrainRegion",
        "RIN",
        "Plate",
        "numReads",
        "ERCCsumLogErr",
        "geneticRace",
        "dropMetrics",
        "dropRegion",
        "dropRace"
    )
pd[which(pd$dropMetrics == TRUE |
             pd$dropRegion == TRUE | pd$dropRace == TRUE), info_cols]

#         BrNum AgeDeath Sex Race PrimaryDx BrainRegion RIN Plate dropMetrics dropRegion dropRace
# R17520 Br1635 52.31000   M CAUC       MDD        sACC 6.6     6       FALSE       TRUE    FALSE
# R17527 Br1675 32.10000   M CAUC       MDD        sACC 6.8     2       FALSE       TRUE    FALSE
# R17538 Br1723 55.70000   F CAUC       MDD    Amygdala 5.3     1        TRUE      FALSE    FALSE
# R17547 Br1754 25.98000   F CAUC       MDD        sACC 5.7     5       FALSE       TRUE    FALSE
# R17579 Br6099 26.66381   M CAUC       MDD        sACC 6.9     7       FALSE       TRUE    FALSE
# R17616 Br1986 25.35000   M CAUC       MDD    Amygdala 6.7     5       FALSE      FALSE     TRUE
# R17794 Br5763 53.89000   M CAUC       MDD        sACC 7.1     1        TRUE      FALSE    FALSE
# R17839 Br5965 30.12183   M CAUC       MDD    Amygdala 5.6     4        TRUE      FALSE    FALSE
# R17880 Br6227 18.62822   M HISP   Control        sACC 7.4     5       FALSE      FALSE     TRUE
# R17894 Br8017 49.29500   M CAUC       MDD        sACC 7.2     3       FALSE       TRUE    FALSE
# R17914 Br8053 47.00616   M CAUC   Control    Amygdala 5.5     1        TRUE      FALSE    FALSE
# R17941 Br8133 62.36813   M CAUC       MDD        sACC 5.7     1       FALSE       TRUE    FALSE
# R17951 Br2333 62.99000   F CAUC   Control        sACC 5.0     4       FALSE       TRUE    FALSE
# R17952 Br5549 67.47000   M CAUC       MDD        sACC 5.6     5       FALSE       TRUE    FALSE
# R17963 Br5526 68.85000   M CAUC       MDD    Amygdala 5.4     4        TRUE      FALSE    FALSE
# R18425 Br6285 64.82000   M   AS   Control    Amygdala 6.8     3       FALSE      FALSE     TRUE
# R18426 Br6285 64.82000   M   AS   Control        sACC 7.2     3       FALSE      FALSE     TRUE
# R18430 Br3863 69.33904   M CAUC   Control        sACC 6.5     7       FALSE       TRUE    FALSE
# R18431 Br3872 59.24162   F HISP   Control    Amygdala 5.8     4       FALSE      FALSE     TRUE
# R18432 Br3872 59.24162   F HISP   Control        sACC 5.9     4       FALSE      FALSE     TRUE
# R18458 Br8313 55.41672   M CAUC   Control        sACC 5.5     4       FALSE       TRUE    FALSE
# R18459 Br6123 30.76249   F   AA   Control    Amygdala 6.6     3       FALSE      FALSE     TRUE
# R18460 Br6123 30.76249   F   AA   Control        sACC 7.1     4       FALSE      FALSE     TRUE
# R18461 Br5146 47.01000   F   AA   Control    Amygdala 6.8     4       FALSE      FALSE     TRUE
# R18462 Br5146 47.01000   F   AA   Control        sACC 7.6     4       FALSE      FALSE     TRUE
# R18824 Br3877 51.85775   M CAUC       MDD    Amygdala 6.4     6       FALSE      FALSE     TRUE
# R18938 Br6293 27.85763   F CAUC       MDD    Amygdala 7.6     5       FALSE       TRUE    FALSE
# R19186 Br1475 43.30000   M CAUC       MDD        sACC 8.9     4       FALSE       TRUE    FALSE
# R13935 Br0922 51.35000   M CAUC   Control    Amygdala 6.2    NA        TRUE      FALSE    FALSE
# R14081 Br5611 63.82000   F CAUC   Bipolar    Amygdala 7.0    NA       FALSE       TRUE    FALSE
# R14110 Br1562 36.15000   M CAUC   Bipolar        sACC 7.1    NA       FALSE       TRUE    FALSE
# R14175 Br1752 39.20000   F CAUC   Bipolar        sACC 7.5    NA       FALSE       TRUE    FALSE
# R14182 Br1558 23.43000   F CAUC   Control        sACC 6.9    NA       FALSE       TRUE    FALSE
# R14210 Br1148 29.97000   M   AA   Control        sACC 7.9    NA       FALSE      FALSE     TRUE
# R14228 Br2333 62.99000   F CAUC   Control        sACC 6.2    NA       FALSE       TRUE    FALSE
# R14235 Br1661 51.15000   M CAUC   Bipolar        sACC 6.2    NA       FALSE       TRUE    FALSE
# R14269 Br5555 21.34000   M CAUC   Control        sACC 6.5    NA       FALSE       TRUE    FALSE
# R14295 Br2071 40.62000   F CAUC   Bipolar        sACC 6.5    NA       FALSE       TRUE    FALSE
# R14306 Br5435 56.51000   M CAUC   Bipolar        sACC 7.9    NA        TRUE      FALSE    FALSE
# R15072 Br5939 55.90965   M CAUC   Bipolar    Amygdala 7.5    NA       FALSE       TRUE    FALSE

#### Summarize ####

pd$dropSum <-
    rowSums(pd[, c("dropMetrics", "dropRegion", "dropRace")])
sum(pd$dropSum > 0) # 40
rownames(pd) <- paste0(pd$RNum, "_", pd$Experiment)

## Save qc_dropping results
qcresults <- pd[, info_cols]

## Update pd
pd <- pd[-which(pd$dropSum > 0),]

message("Remaining Samples: n", nrow(pd))
# 1100

table(pd$BrainRegion, pd$PrimaryDx)
table(pd$Sex, pd$PrimaryDx)
table(pd$Race, pd$PrimaryDx)
table(table(pd$BrNum))
summary(pd$Age)


# > table(pd$BrainRegion, pd$PrimaryDx)
#
# Bipolar Control MDD
# Amygdala     123     188 234
# sACC         123     202 230
# > table(pd$Sex, pd$PrimaryDx)
#
# Bipolar Control MDD
# F      96      77 156
# M     150     313 308
# > table(pd$Race, pd$PrimaryDx)
#
# Bipolar Control MDD
# CAUC     246     389 464
# HISP       0       1   0
# > table(table(pd$BrNum))
#
# 1   2   3
# 98 492   6
# > summary(pd$Age)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 17.37   34.53   47.09   46.55   55.87   95.27

rownames(qcresults) <- paste0(qcresults$RNum,"_",qcresults$Experiment)
write.csv(qcresults, file = "qc_dropping_results.csv")


#### PCA - Combined data ####
rse_gene <- rse_gene[, rse_gene$RNum %in% pd$RNum]
gRpkmBP <- recount::getRPKM(rse_gene, "Length")
## filter low expressing genes
gRpkmBP <- gRpkmBP[which(rowMeans(gRpkmBP) > 0.2),]
yExprsComb <- log2(gRpkmBP + 1)

# combined phenodata
pd$group <- paste0(pd$BrainRegion, "_", pd$PrimaryDx)
pca1 <- prcomp(t(yExprsComb))
pcaVars1 <- getPcaVars(pca1)

pdf("pdfs/pca_log2Rpkm_PC1_2_combined_datasets.pdf",
    h = 6,
    w = 6)
par(
    mar = c(5, 6, 4, 2),
    cex.axis = 1.5,
    cex.lab = 1.5,
    cex.main = 1.5
)
palette(brewer.pal(8, "Spectral"))
plot(
    pca1$x,
    pch = 21,
    bg = (as.numeric(factor(pd$Experiment)) * 4) - 1,
    cex = 2,
    main = "Gene PCs",
    xlab = paste0("PC1: ", pcaVars1[1], "% Var Expl"),
    ylab = paste0("PC2: ", pcaVars1[2], "% Var Expl")
)
legend(
    "topleft",
    paste0(levels(factor(pd$Experiment))),
    pch = 15,
    col = c(3, 7),
    cex = .9
)
plot(
    pca1$x,
    pch = 21,
    bg = as.numeric(factor(pd$PrimaryDx)) * 2,
    cex = 2,
    main = "Gene PCs",
    xlab = paste0("PC1: ", pcaVars1[1], "% Var Expl"),
    ylab = paste0("PC2: ", pcaVars1[2], "% Var Expl")
)
legend(
    "topleft",
    paste0(levels(factor(pd$PrimaryDx))),
    pch = 15,
    col = c(2, 4, 6),
    cex = .9
)
plot(
    pca1$x,
    pch = 21,
    bg = as.numeric(factor(pd$BrainRegion)) * 3,
    cex = 2,
    main = "Gene PCs",
    xlab = paste0("PC1: ", pcaVars1[1], "% Var Expl"),
    ylab = paste0("PC2: ", pcaVars1[2], "% Var Expl")
)
legend(
    "topleft",
    paste0(levels(factor(pd$BrainRegion))),
    pch = 15,
    col = c(3, 6),
    cex = .9
)
plot(
    pca1$x,
    pch = 21,
    bg = as.numeric(factor(pd$group)),
    cex = 2,
    main = "Gene PCs",
    xlab = paste0("PC1: ", pcaVars1[1], "% Var Expl"),
    ylab = paste0("PC2: ", pcaVars1[2], "% Var Expl")
)
legend(
    "topleft",
    paste0(levels(factor(pd$group))),
    pch = 15,
    col = 1:6,
    cex = .9
)
dev.off()

#### PCA - Filter for MDD samples ####
rse_gene <- rse_gene[, rse_gene$Experiment == "psychENCODE_MDD"]
pd_mdd <- pd[pd$Experiment == "psychENCODE_MDD",]
gRpkm <- recount::getRPKM(rse_gene, "Length")

## filter low expressing genes
gRpkm <- gRpkm[which(rowMeans(gRpkm) > 0.2),]
yExprs <- log2(gRpkm + 1)

pca1 <- prcomp(t(yExprs))
pcaVars1 <- getPcaVars(pca1)

pdf("pdfs/pca_log2Rpkm_PC1_2.pdf", h = 6, w = 6)
par(
    mar = c(5, 6, 4, 2),
    cex.axis = 1.5,
    cex.lab = 1.5,
    cex.main = 1.5
)
palette(brewer.pal(8, "Spectral"))
plot(
    pca1$x,
    pch = 21,
    bg = as.numeric(factor(pd_mdd$PrimaryDx)) * 3,
    cex = 2,
    main = "Gene PCs",
    xlab = paste0("PC1: ", pcaVars1[1], "% Var Expl"),
    ylab = paste0("PC2: ", pcaVars1[2], "% Var Expl")
)
legend(
    "topleft",
    paste0(levels(factor(pd_mdd$PrimaryDx))),
    pch = 15,
    col = c(3, 6),
    cex = .9
)
plot(
    pca1$x,
    pch = 21,
    bg = as.numeric(factor(pd_mdd$Sex)) * 2,
    cex = 2,
    main = "Gene PCs",
    xlab = paste0("PC1: ", pcaVars1[1], "% Var Expl"),
    ylab = paste0("PC2: ", pcaVars1[2], "% Var Expl")
)
legend(
    "topleft",
    paste0(levels(factor(pd_mdd$Sex))),
    pch = 15,
    col = c(2, 4),
    cex = .9
)
plot(
    pca1$x,
    pch = 21,
    bg = as.numeric(factor(pd_mdd$BrainRegion)) * 3,
    cex = 2,
    main = "Gene PCs",
    xlab = paste0("PC1: ", pcaVars1[1], "% Var Expl"),
    ylab = paste0("PC2: ", pcaVars1[2], "% Var Expl")
)
legend(
    "topleft",
    paste0(levels(factor(pd_mdd$BrainRegion))),
    pch = 15,
    col = c(3, 6),
    cex = .9
)
plot(
    pca1$x,
    pch = 21,
    bg = as.numeric(factor(pd_mdd$Race)) * 2,
    cex = 2,
    main = "Gene PCs",
    xlab = paste0("PC1: ", pcaVars1[1], "% Var Expl"),
    ylab = paste0("PC2: ", pcaVars1[2], "% Var Expl")
)
legend(
    "topleft",
    paste0(levels(factor(pd_mdd$Race))),
    pch = 15,
    col = c(2, 4),
    cex = .9
)
plot(
    pca1$x,
    pch = 21,
    bg = as.numeric(factor(pd_mdd$Plate)),
    cex = 2,
    main = "Gene PCs",
    xlab = paste0("PC1: ", pcaVars1[1], "% Var Expl"),
    ylab = paste0("PC2: ", pcaVars1[2], "% Var Expl")
)
legend(
    "topleft",
    paste0("Plate ", levels(factor(pd_mdd$Plate))),
    pch = 15,
    col = 1:8,
    cex = .9
)
dev.off()


# sgejobs::job_single('qc_checks', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript qc_checks.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
