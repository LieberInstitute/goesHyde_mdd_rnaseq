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
library(here)

## For styling this script
# styler::style_file("qc_checks.R", transformers = biocthis::bioc_style())

## load phenotype and alignment data
load(here("data","rse_gene_raw_GoesZandi.rda"),
     verbose = TRUE
)
pd <- colData(rse_gene) %>% as.data.frame()
table(pd$overlap)
# FALSE  TRUE 
# 1123    20 

metricCols <-c("RNum", "Experiment", "numReads", "Plate", "ERCCsumLogErr")
pd_mdd <- read.csv(here("preprocessed_data","read_and_alignment_metrics_goesHyde_MDD.csv"),
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
rownames(pd) <- paste0(pd$RNum,"_",pd$Experiment)

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
abline(v = 1e-3, lty = 2)
legend("bottomright", "c)", bty = "n", cex = 2)
plot(pd$RIN, log10(pd$numReads), pch = 21, bg = "grey")
abline(h = log10(10^7.25), lty = 2)
legend("bottomright", "d)", bty = "n", cex = 2)
dev.off()


####################################

## drop samples
pd$dropMetrics <- FALSE
pd$dropMetrics[pd$overallMapRate < 0.5 |
                   pd$totalAssignedGene < .3 | 
                   pd$numReads < (10^7.25) |
                   pd$rRNA_rate > 1e-3] <- TRUE

table(pd$dropMetrics)
# FALSE  TRUE
# 1133     10

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
abline(v = 1e-3, lty = 2)
plot(
    pd$RIN,
    log10(pd$numReads),
    pch = 21,
    bg = pd$dropMetrics + 1,
    cex = pd$dropMetrics + 1,
    ylim = c(log10(rL), log10(rH))
)
abline(h = log10((10^7.25)), lty = 2)
dev.off()


#### ERCC vs. Metrics ####
pdf("pdfs/ERCC_check_postdrop.pdf", h = 10, w = 10)
par(
    mfcol = c(2, 2),
    mar = c(5, 6, 2, 2),
    cex.axis = 1.8,
    cex.lab = 1.8
)
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
abline(h = log10(10^7.25), lty = 2)
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
abline(h = 1e-3, lty = 2)
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
stopifnot(identical(names(guess), rownames(pd)))
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


# BrNum      Experiment AgeDeath Sex Race PrimaryDx BrainRegion dropMetrics       guess
# R14179_psychENCODE_MDD Br1469 psychENCODE_MDD 28.57000   M CAUC   Control        sACC       FALSE -0.01452763
# R17520_psychENCODE_MDD Br1635 psychENCODE_MDD 52.31000   M CAUC       MDD        sACC       FALSE  0.20640366
# R17527_psychENCODE_MDD Br1675 psychENCODE_MDD 32.10000   M CAUC       MDD        sACC       FALSE -0.03900185
# R17547_psychENCODE_MDD Br1754 psychENCODE_MDD 25.98000   F CAUC       MDD        sACC       FALSE  0.02626940
# R17579_psychENCODE_MDD Br6099 psychENCODE_MDD 26.66381   M CAUC       MDD        sACC       FALSE  0.12403490
# R17894_psychENCODE_MDD Br8017 psychENCODE_MDD 49.29500   M CAUC       MDD        sACC       FALSE  0.25591506
# R17941_psychENCODE_MDD Br8133 psychENCODE_MDD 62.36813   M CAUC       MDD        sACC       FALSE -0.09759443
# R17952_psychENCODE_MDD Br5549 psychENCODE_MDD 67.47000   M CAUC       MDD        sACC       FALSE  0.01866621
# R18430_psychENCODE_MDD Br3863 psychENCODE_MDD 69.33904   M CAUC   Control        sACC       FALSE -0.12428720
# R18458_psychENCODE_MDD Br8313 psychENCODE_MDD 55.41672   M CAUC   Control        sACC       FALSE  0.14767064
# R19186_psychENCODE_MDD Br1475 psychENCODE_MDD 43.30000   M CAUC       MDD        sACC       FALSE  0.04996062
# R14110_psychENCODE_BP  Br1562  psychENCODE_BP 36.15000   M CAUC   Bipolar        sACC       FALSE  0.27638628
# R14175_psychENCODE_BP  Br1752  psychENCODE_BP 39.20000   F CAUC   Bipolar        sACC       FALSE  0.06713448
# R14179_psychENCODE_BP  Br1469  psychENCODE_BP 28.57000   M CAUC   Control        sACC       FALSE -0.03493316
# R14182_psychENCODE_BP  Br1558  psychENCODE_BP 23.43000   F CAUC   Control        sACC       FALSE  0.26839632
# R14228_psychENCODE_BP  Br2333  psychENCODE_BP 62.99000   F CAUC   Control        sACC       FALSE  0.03937358
# R14235_psychENCODE_BP  Br1661  psychENCODE_BP 51.15000   M CAUC   Bipolar        sACC       FALSE  0.27786712
# R14269_psychENCODE_BP  Br5555  psychENCODE_BP 21.34000   M CAUC   Control        sACC       FALSE  0.13650520
# R14295_psychENCODE_BP  Br2071  psychENCODE_BP 40.62000   F CAUC   Bipolar        sACC       FALSE  0.14280649
# > pd[which(pd$BrainRegion == "Amygdala" & pd$guess > 0.9), info_cols]
# BrNum      Experiment AgeDeath Sex Race PrimaryDx BrainRegion dropMetrics    guess
# R17496_psychENCODE_MDD Br1469 psychENCODE_MDD 28.57000   M CAUC   Control    Amygdala       FALSE 1.073050
# R14081_psychENCODE_BP  Br5611  psychENCODE_BP 63.82000   F CAUC   Bipolar    Amygdala       FALSE 1.014568
# R15072_psychENCODE_BP  Br5939  psychENCODE_BP 55.90965   M CAUC   Bipolar    Amygdala       FALSE 1.055168

## Br1469 labels got switched
pd[which(pd$RNum == "R14179"), "BrainRegion"] <- "Amygdala"
pd[which(pd$RNum == "R17496"), "BrainRegion"] <- "sACC"

## Drop the rest
dropInd <- which((pd$BrainRegion == "sACC" & pd$guess < .3) |
                     (pd$BrainRegion == "Amygdala" & pd$guess > 0.9))
pd$dropRegion <- FALSE
pd$dropRegion[dropInd] <- TRUE

table(pd$dropRegion)
# FALSE  TRUE
# 1119    21
## With overlap controls included
# FALSE  TRUE 
# 1126    19

######################################################## .
gia <- read.csv(here("synapse","genotypeInferredAncestry.csv"))
pd$geneticRace <- pd$Race
pd$geneticRace[match(gia$BrNum, pd$BrNum)] <-
    gia$genotypeInferredAncestry

pd$dropRace <- pd$geneticRace != "CAUC"

table(pd$Race, pd$geneticRace)
# AA ambiguous   AS CAUC HISP
# AA      5         0    0    0    0
# AS      0         0    2    0    0
# CAUC    0         2    0 1130    0
# HISP    0         0    0    1    3

table(pd$dropRace)
# FALSE  TRUE 
# 1131    12 


######################################################## .
temp_info_cols <-
    c(  "BrNum",
        # "AgeDeath",
        "Sex",
        "Race",
        "PrimaryDx",
        "BrainRegion",
        "RIN",
        # "Plate",
        "numReads",
        # "ERCCsumLogErr",
        "geneticRace",
        "dropMetrics",
        "dropRegion",
        "dropRace"
    )
pd[which(pd$dropMetrics == TRUE |
             pd$dropRegion == TRUE | pd$dropRace == TRUE), temp_info_cols]

# BrNum Sex Race PrimaryDx BrainRegion RIN  numReads geneticRace dropMetrics dropRegion dropRace
# R17514_psychENCODE_MDD Br1582   F CAUC       MDD    Amygdala 6.7  15252860        CAUC        TRUE      FALSE    FALSE
# R17520_psychENCODE_MDD Br1635   M CAUC       MDD        sACC 6.6  72006392        CAUC       FALSE       TRUE    FALSE
# R17527_psychENCODE_MDD Br1675   M CAUC       MDD        sACC 6.8  95710478        CAUC       FALSE       TRUE    FALSE
# R17538_psychENCODE_MDD Br1723   F CAUC       MDD    Amygdala 5.3  49607375        CAUC        TRUE      FALSE    FALSE
# R17547_psychENCODE_MDD Br1754   F CAUC       MDD        sACC 5.7 164708040        CAUC       FALSE       TRUE    FALSE
# R17579_psychENCODE_MDD Br6099   M CAUC       MDD        sACC 6.9  84212528        CAUC       FALSE       TRUE    FALSE
# R17616_psychENCODE_MDD Br1986   M CAUC       MDD    Amygdala 6.7 325417972   ambiguous       FALSE      FALSE     TRUE
# R17794_psychENCODE_MDD Br5763   M CAUC       MDD        sACC 7.1 132887932        CAUC        TRUE      FALSE    FALSE
# R17839_psychENCODE_MDD Br5965   M CAUC       MDD    Amygdala 5.6 111654297        CAUC        TRUE      FALSE    FALSE
# R17880_psychENCODE_MDD Br6227   M HISP   Control        sACC 7.4 110688364        HISP       FALSE      FALSE     TRUE
# R17894_psychENCODE_MDD Br8017   M CAUC       MDD        sACC 7.2 159142448        CAUC       FALSE       TRUE    FALSE
# R17914_psychENCODE_MDD Br8053   M CAUC   Control    Amygdala 5.5 216752444        CAUC        TRUE      FALSE    FALSE
# R17941_psychENCODE_MDD Br8133   M CAUC       MDD        sACC 5.7  88944792        CAUC       FALSE       TRUE    FALSE
# R17950_psychENCODE_MDD Br1698   M CAUC       MDD        sACC 6.6  10284026        CAUC        TRUE      FALSE    FALSE
# R17952_psychENCODE_MDD Br5549   M CAUC       MDD        sACC 5.6 125112170        CAUC       FALSE       TRUE    FALSE
# R17963_psychENCODE_MDD Br5526   M CAUC       MDD    Amygdala 5.4 192819932        CAUC        TRUE      FALSE    FALSE
# R18425_psychENCODE_MDD Br6285   M   AS   Control    Amygdala 6.8  63210188          AS       FALSE      FALSE     TRUE
# R18426_psychENCODE_MDD Br6285   M   AS   Control        sACC 7.2 191881132          AS       FALSE      FALSE     TRUE
# R18430_psychENCODE_MDD Br3863   M CAUC   Control        sACC 6.5  94282744        CAUC       FALSE       TRUE    FALSE
# R18431_psychENCODE_MDD Br3872   F HISP   Control    Amygdala 5.8 208135920        HISP       FALSE      FALSE     TRUE
# R18432_psychENCODE_MDD Br3872   F HISP   Control        sACC 5.9 135427576        HISP       FALSE      FALSE     TRUE
# R18458_psychENCODE_MDD Br8313   M CAUC   Control        sACC 5.5 226538640        CAUC       FALSE       TRUE    FALSE
# R18459_psychENCODE_MDD Br6123   F   AA   Control    Amygdala 6.6 102404724          AA       FALSE      FALSE     TRUE
# R18460_psychENCODE_MDD Br6123   F   AA   Control        sACC 7.1 223245052          AA       FALSE      FALSE     TRUE
# R18461_psychENCODE_MDD Br5146   F   AA   Control    Amygdala 6.8  89905868          AA       FALSE      FALSE     TRUE
# R18462_psychENCODE_MDD Br5146   F   AA   Control        sACC 7.6 179353610          AA       FALSE      FALSE     TRUE
# R18824_psychENCODE_MDD Br3877   M CAUC       MDD    Amygdala 6.4 125328900   ambiguous       FALSE      FALSE     TRUE
# R19186_psychENCODE_MDD Br1475   M CAUC       MDD        sACC 8.9 313300392        CAUC       FALSE       TRUE    FALSE
# R13935_psychENCODE_BP  Br0922   M CAUC   Control    Amygdala 6.2 116357154        CAUC        TRUE      FALSE    FALSE
# R14003_psychENCODE_BP  Br2333   F CAUC   Control    Amygdala 6.4 139167058        CAUC        TRUE      FALSE    FALSE
# R14081_psychENCODE_BP  Br5611   F CAUC   Bipolar    Amygdala 7.0 133289152        CAUC       FALSE       TRUE    FALSE
# R14110_psychENCODE_BP  Br1562   M CAUC   Bipolar        sACC 7.1 137413506        CAUC       FALSE       TRUE    FALSE
# R14175_psychENCODE_BP  Br1752   F CAUC   Bipolar        sACC 7.5 139218132        CAUC       FALSE       TRUE    FALSE
# R14182_psychENCODE_BP  Br1558   F CAUC   Control        sACC 6.9 137971714        CAUC       FALSE       TRUE    FALSE
# R14210_psychENCODE_BP  Br1148   M   AA   Control        sACC 7.9 144816192          AA       FALSE      FALSE     TRUE
# R14228_psychENCODE_BP  Br2333   F CAUC   Control        sACC 6.2 130240626        CAUC       FALSE       TRUE    FALSE
# R14235_psychENCODE_BP  Br1661   M CAUC   Bipolar        sACC 6.2 123042494        CAUC       FALSE       TRUE    FALSE
# R14269_psychENCODE_BP  Br5555   M CAUC   Control        sACC 6.5 119019314        CAUC       FALSE       TRUE    FALSE
# R14295_psychENCODE_BP  Br2071   F CAUC   Bipolar        sACC 6.5 134903940        CAUC       FALSE       TRUE    FALSE
# R14306_psychENCODE_BP  Br5435   M CAUC   Bipolar        sACC 7.9 186118990        CAUC        TRUE      FALSE    FALSE
# R15072_psychENCODE_BP  Br5939   M CAUC   Bipolar    Amygdala 7.5 153561278        CAUC       FALSE       TRUE    FALSE

#### Summarize ####
pd$dropSum <-
    rowSums(pd[, c("dropMetrics", "dropRegion", "dropRace")])
message("Number of samples to drop: ",sum(pd$dropSum > 0))

## Save qc_dropping results
qcresults <- pd

## Update pd
pd <- pd[-which(pd$dropSum > 0 & !(pd$overlap & pd$Experiment == "psychENCODE_BP")),]

message("Remaining Samples: n", nrow(pd))
# 1102

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
