####

### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)
library(RColorBrewer)
library(sessioninfo)
library(here)

## load
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_exon.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_jxn.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_tx.Rdata"), verbose = TRUE)


pd <- colData(rse_gene)

## load SNP data
load(here("eqtl", "genomewide", "rdas", "overlappingSNPs.rda"), verbose = TRUE) # snpMapKeep
load(here("genotype_data", "goesHyde_bipolarMdd_Genotypes_n593.rda"), verbose = TRUE) # snp, snpMap & mds
snpMap$pos_hg19 <- paste0(snpMap$CHR, ":", snpMap$POS)

#  In Bipolar code rownames are SNP, need to assign here
# > identical(rownames(snpMap), snpMap$SNP)
# [1] TRUE
rownames(snp) <- rownames(snpMap) <- snpMap$SNP
snpInd <- which(rownames(snpMap) %in% rownames(snpMapKeep) & !is.na(snpMap$pos_hg38))
snpMap <- snpMap[snpInd, ]
snp <- snp[snpInd, ]

#####################
# filter brain region
# make mds and snp dimensions equal to N
# (repeat rows or columns for BrNum replicates)
table(pd$BrNum %in% colnames(snp))
table(pd$BrNum %in% rownames(mds))

has_genotype <- pd$BrNum %in% rownames(mds) # snp and mds are missing BR1843
rse_gene <- rse_gene[, has_genotype]
rse_exon <- rse_exon[, has_genotype]
rse_jxn <- rse_jxn[, has_genotype]
rse_tx <- rse_tx[, has_genotype]

pd <- pd[has_genotype, ]
mds <- mds[pd$BrNum, ]
snp <- snp[, pd$BrNum]
rownames(mds) <- colnames(snp) <- pd$RNum

snpMap$maf <- rowSums(snp, na.rm = TRUE) / (2 * rowSums(!is.na(snp)))


################
## load table ##
###############
load("mergedEqtl_output_sacc_genomewide_4features_FDR01.rda", verbose = TRUE)
sacc <- allEqtlFDR01

pd$PrimaryDx <- factor(pd$PrimaryDx,
    levels = c("Control", "Bipolar", "MDD")
)
mod <- model.matrix(~ PrimaryDx + Sex + as.matrix(mds[, 1:5]), data = pd)

#####################
## calculate rpkm ##
####################
message(Sys.time(), " Get rpkm values")
geneRpkm <- recount::getRPKM(rse_gene, "Length")
exonRpkm <- recount::getRPKM(rse_exon, "Length")
rowData(rse_jxn)$Length <- 100
jxnRp10m <- recount::getRPKM(rse_jxn, "Length")
txTpm <- assays(rse_tx)$tpm

## residualize expression
gExprs <- log2(geneRpkm + 1)
gExprs <- cleaningY(gExprs, mod, P = 1)

eExprs <- log2(exonRpkm + 1)
eExprs <- cleaningY(eExprs, mod, P = 1)

jExprs <- log2(jxnRp10m + 1)
jExprs <- cleaningY(jExprs, mod, P = 1)

tExprs <- log2(txTpm + 1)
tExprs <- cleaningY(tExprs, mod, P = 1)


exprsAdj <- rbind(gExprs, eExprs, jExprs, tExprs)
sacc$Symbol <- as.character(sacc$Symbol)

saccG <- sacc[which(sacc$Type == "Gene"), ]
saccE <- sacc[which(sacc$Type == "Exon"), ]
saccJ <- sacc[which(sacc$Type == "Jxn"), ]
saccT <- sacc[which(sacc$Type == "Tx"), ]

##########
## plot ##
#########
message(Sys.time(), " Box plots")
pdf("sacc_top_eqtl_adj.pdf", h = 6, w = 10)
par(mfrow = c(2, 3), cex.main = 1.2, cex.lab = 1.2)
palette(brewer.pal(8, "Spectral"))

for (i in 1:12) {
    symi <- saccG[i, "Symbol"]
    symi[is.na(symi)] <- ""
    snpi <- saccG[i, "snps"]
    feati <- saccG[i, "gene"]
    p_i <- signif(saccG[i, "pvalue"], 3)
    typei <- saccG[i, "Type"]

    boxplot(exprsAdj[feati, ] ~ unlist(snp[snpi, ]),
        xlab = snpi, ylab = "Residualized Expression",
        main = paste0(symi, "\n", feati, " (", typei, ")"),
        ylim = c(range(exprsAdj[feati, ])), outline = FALSE
    )
    points(exprsAdj[feati, ] ~ jitter(unlist(snp[snpi, ]) + 1),
        pch = 21,
        bg = as.numeric(unlist(snp[snpi, ])) + 2, cex = 1.5
    )
    legend("top", paste0("p=", p_i))
}
for (i in 1:12) {
    symi <- saccE[i, "Symbol"]
    symi[is.na(symi)] <- ""
    snpi <- saccE[i, "snps"]
    feati <- saccE[i, "gene"]
    p_i <- signif(saccE[i, "pvalue"], 3)
    typei <- saccE[i, "Type"]

    boxplot(exprsAdj[feati, ] ~ unlist(snp[snpi, ]),
        xlab = snpi, ylab = "Residualized Expression",
        main = paste0(symi, "\n", feati, " (", typei, ")"),
        ylim = c(range(exprsAdj[feati, ])), outline = FALSE
    )
    points(exprsAdj[feati, ] ~ jitter(unlist(snp[snpi, ]) + 1),
        pch = 21,
        bg = as.numeric(unlist(snp[snpi, ])) + 2, cex = 1.5
    )
    legend("top", paste0("p=", p_i))
}
for (i in 1:12) {
    symi <- saccJ[i, "Symbol"]
    symi[is.na(symi)] <- ""
    snpi <- saccJ[i, "snps"]
    feati <- saccJ[i, "gene"]
    p_i <- signif(saccJ[i, "pvalue"], 3)
    typei <- saccJ[i, "Type"]

    boxplot(exprsAdj[feati, ] ~ unlist(snp[snpi, ]),
        xlab = snpi, ylab = "Residualized Expression",
        main = paste0(symi, "\n", feati, " (", typei, ")"),
        ylim = c(range(exprsAdj[feati, ])), outline = FALSE
    )
    points(exprsAdj[feati, ] ~ jitter(unlist(snp[snpi, ]) + 1),
        pch = 21,
        bg = as.numeric(unlist(snp[snpi, ])) + 2, cex = 1.5
    )
    legend("top", paste0("p=", p_i))
}
for (i in 1:12) {
    symi <- saccT[i, "Symbol"]
    symi[is.na(symi)] <- ""
    snpi <- saccT[i, "snps"]
    feati <- saccT[i, "gene"]
    p_i <- signif(saccT[i, "pvalue"], 3)
    typei <- saccT[i, "Type"]

    boxplot(exprsAdj[feati, ] ~ unlist(snp[snpi, ]),
        xlab = snpi, ylab = "Residualized Expression",
        main = paste0(symi, "\n", feati, " (", typei, ")"),
        ylim = c(range(exprsAdj[feati, ])), outline = FALSE
    )
    points(exprsAdj[feati, ] ~ jitter(unlist(snp[snpi, ]) + 1),
        pch = 21,
        bg = as.numeric(unlist(snp[snpi, ])) + 2, cex = 1.5
    )
    legend("top", paste0("p=", p_i))
}
dev.off()

##################################
## make CSV of top 1000 of each ##
#################################
message(Sys.time(), " Create csv")

sacc <- rbind(saccG[1:1000, ], saccE[1:1000, ], saccJ[1:1000, ], saccT[1:1000, ])
remove(saccG, saccE, saccJ, saccT)
sacc$EnsemblGeneID <- ss(sacc$EnsemblGeneID, "\\.")

## snpMap
# load(here("genotype_data","goesHyde_bipolarMdd_Genotypes_n593.rda"), verbose = TRUE) # do we need to reload this?
# snpMap$hg19POS = paste0(snpMap$CHR,":",snpMap$POS)
# # snpMap = snpMap[which(rownames(snpMap) %in% c(sacc$snps,sacc$snps,dlp$snps) ),c("SNP","chr_hg38","pos_hg38","hg19POS")]

## featMap
load(here("data", "rse_gene_GoesZandi_n1093.rda"), verbose = TRUE)
load(here("data", "rse_exon_GoesZandi_n1093.rda"), verbose = TRUE)
load(here("data", "rse_jxn_GoesZandi_n1093.rda"), verbose = TRUE)
load(here("data", "rse_tx_GoesZandi_n1093.rda"), verbose = TRUE)
gMap <- as.data.frame(rowRanges(rse_gene))[, c("seqnames", "start", "end", "strand", "Class")]
eMap <- as.data.frame(rowRanges(rse_exon))[, c("seqnames", "start", "end", "strand", "Class")]
jMap <- as.data.frame(rowRanges(rse_jxn))[, c("seqnames", "start", "end", "strand", "Class")]
txMap <- as.data.frame(rowRanges(rse_tx))[, c("seqnames", "start", "end", "strand", "source")]
txMap$source <- "InGen"
# rm(rse_gene, rse_exon, rse_jxn, rse_tx)
colnames(gMap) <- colnames(eMap) <- colnames(jMap) <- colnames(txMap) <-
    c("feat_chr", "feat_start", "feat_end", "strand", "Class")
featMap <- rbind(rbind(rbind(gMap, eMap), jMap), txMap)
featMap$Type <- c(rep("Gene", nrow(gMap)), rep("Exon", nrow(eMap)), rep("Jxn", nrow(jMap)), rep("Tx", nrow(txMap)))

geneMap <- as.data.frame(rowRanges(rse_gene))[, c("gencodeID", "Symbol", "ensemblID", "gene_type")]

## put together
snpMap_temp <- snpMap[sacc$snps, ]
featMap_temp <- featMap[sacc$gene, ]
geneMap_temp <- geneMap[match(sacc$EnsemblGeneID, geneMap$ensemblID), ]
sacc2 <- cbind(cbind(cbind(snpMap_temp, featMap_temp), geneMap_temp), sacc)
# sacc2 = sacc2[,-which(colnames(sacc2)=="gencodeTx")]

sacc3 <- sacc2[, c(2, 12:14, 26, 20, 15:19, 22:24, 27:30)]
write.csv(sacc3, file = "genomewide_snps_sacc_eqtls_top1000.csv")


# sgejobs::job_single("make_boxplots_and_csv_sacc", memory = "150G",create_shell = TRUE, command = "Rscript make_boxplots_and_csv_sacc.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# [1] "2020-06-08 18:08:50 EDT"
#     user   system  elapsed
# 3032.423  653.619 3692.820
# ??? Session info ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#  setting  value
#  version  R version 4.0.0 Patched (2020-05-13 r78451)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2020-06-08
#
# ??? Packages ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#  package              * version  date       lib source
#  acepack                1.4.1    2016-10-29 [2] CRAN (R 4.0.0)
#  annotate               1.66.0   2020-04-27 [2] Bioconductor
#  AnnotationDbi          1.50.0   2020-04-27 [2] Bioconductor
#  askpass                1.1      2019-01-13 [2] CRAN (R 4.0.0)
#  assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.0.0)
#  backports              1.1.7    2020-05-13 [2] CRAN (R 4.0.0)
#  base64enc              0.1-3    2015-07-28 [2] CRAN (R 4.0.0)
#  Biobase              * 2.48.0   2020-04-27 [2] Bioconductor
#  BiocFileCache          1.12.0   2020-04-27 [2] Bioconductor
#  BiocGenerics         * 0.34.0   2020-04-27 [2] Bioconductor
#  BiocParallel         * 1.22.0   2020-04-27 [2] Bioconductor
#  biomaRt                2.44.0   2020-04-27 [2] Bioconductor
#  Biostrings             2.56.0   2020-04-27 [2] Bioconductor
#  bit                    1.1-15.2 2020-02-10 [2] CRAN (R 4.0.0)
#  bit64                  0.9-7    2017-05-08 [2] CRAN (R 4.0.0)
#  bitops                 1.0-6    2013-08-17 [2] CRAN (R 4.0.0)
#  blob                   1.2.1    2020-01-20 [2] CRAN (R 4.0.0)
#  BSgenome               1.56.0   2020-04-27 [2] Bioconductor
#  bumphunter             1.30.0   2020-04-27 [2] Bioconductor
#  checkmate              2.0.0    2020-02-06 [2] CRAN (R 4.0.0)
#  cli                    2.0.2    2020-02-28 [1] CRAN (R 4.0.0)
#  cluster                2.1.0    2019-06-19 [3] CRAN (R 4.0.0)
#  codetools              0.2-16   2018-12-24 [3] CRAN (R 4.0.0)
#  colorspace             1.4-1    2019-03-18 [2] CRAN (R 4.0.0)
#  crayon                 1.3.4    2017-09-16 [2] CRAN (R 4.0.0)
#  curl                   4.3      2019-12-02 [2] CRAN (R 4.0.0)
#  data.table             1.12.8   2019-12-09 [2] CRAN (R 4.0.0)
#  DBI                    1.1.0    2019-12-15 [2] CRAN (R 4.0.0)
#  dbplyr                 1.4.3    2020-04-19 [2] CRAN (R 4.0.0)
#  DelayedArray         * 0.14.0   2020-04-27 [2] Bioconductor
#  derfinder              1.22.0   2020-04-27 [2] Bioconductor
#  derfinderHelper        1.22.0   2020-04-27 [2] Bioconductor
#  digest                 0.6.25   2020-02-23 [1] CRAN (R 4.0.0)
#  doRNG                  1.8.2    2020-01-27 [2] CRAN (R 4.0.0)
#  downloader             0.4      2015-07-09 [2] CRAN (R 4.0.0)
#  dplyr                  0.8.5    2020-03-07 [1] CRAN (R 4.0.0)
#  edgeR                  3.30.0   2020-04-27 [2] Bioconductor
#  ellipsis               0.3.1    2020-05-15 [2] CRAN (R 4.0.0)
#  fansi                  0.4.1    2020-01-08 [2] CRAN (R 4.0.0)
#  foreach                1.5.0    2020-03-30 [2] CRAN (R 4.0.0)
#  foreign                0.8-79   2020-04-26 [3] CRAN (R 4.0.0)
#  Formula                1.2-3    2018-05-03 [2] CRAN (R 4.0.0)
#  genefilter           * 1.70.0   2020-04-27 [2] Bioconductor
#  GenomeInfoDb         * 1.24.0   2020-04-27 [2] Bioconductor
#  GenomeInfoDbData       1.2.3    2020-05-18 [2] Bioconductor
#  GenomicAlignments      1.24.0   2020-04-27 [2] Bioconductor
#  GenomicFeatures        1.40.0   2020-04-27 [2] Bioconductor
#  GenomicFiles           1.24.0   2020-04-27 [2] Bioconductor
#  GenomicRanges        * 1.40.0   2020-04-27 [2] Bioconductor
#  GEOquery               2.56.0   2020-04-27 [2] Bioconductor
#  ggplot2                3.3.0    2020-03-05 [2] CRAN (R 4.0.0)
#  glue                   1.4.1    2020-05-13 [1] CRAN (R 4.0.0)
#  googledrive            1.0.1    2020-05-05 [1] CRAN (R 4.0.0)
#  gridExtra              2.3      2017-09-09 [2] CRAN (R 4.0.0)
#  gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.0.0)
#  here                 * 0.1      2017-05-28 [1] CRAN (R 4.0.0)
#  Hmisc                  4.4-0    2020-03-23 [2] CRAN (R 4.0.0)
#  hms                    0.5.3    2020-01-08 [2] CRAN (R 4.0.0)
#  htmlTable              1.13.3   2019-12-04 [2] CRAN (R 4.0.0)
#  htmltools              0.4.0    2019-10-04 [2] CRAN (R 4.0.0)
#  htmlwidgets            1.5.1    2019-10-08 [2] CRAN (R 4.0.0)
#  httr                   1.4.1    2019-08-05 [2] CRAN (R 4.0.0)
#  IRanges              * 2.22.1   2020-04-28 [2] Bioconductor
#  iterators              1.0.12   2019-07-26 [2] CRAN (R 4.0.0)
#  jaffelab             * 0.99.30  2020-05-21 [1] Github (LieberInstitute/jaffelab@42637ff)
#  jpeg                   0.1-8.1  2019-10-24 [2] CRAN (R 4.0.0)
#  jsonlite               1.6.1    2020-02-02 [2] CRAN (R 4.0.0)
#  knitr                  1.28     2020-02-06 [2] CRAN (R 4.0.0)
#  lattice                0.20-41  2020-04-02 [3] CRAN (R 4.0.0)
#  latticeExtra           0.6-29   2019-12-19 [2] CRAN (R 4.0.0)
#  lifecycle              0.2.0    2020-03-06 [1] CRAN (R 4.0.0)
#  limma                  3.44.1   2020-04-28 [2] Bioconductor
#  locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.0.0)
#  magrittr               1.5      2014-11-22 [2] CRAN (R 4.0.0)
#  Matrix                 1.2-18   2019-11-27 [3] CRAN (R 4.0.0)
#  MatrixEQTL           * 2.3      2019-12-22 [1] CRAN (R 4.0.0)
#  matrixStats          * 0.56.0   2020-03-13 [2] CRAN (R 4.0.0)
#  memoise                1.1.0    2017-04-21 [2] CRAN (R 4.0.0)
#  mgcv                 * 1.8-31   2019-11-09 [3] CRAN (R 4.0.0)
#  munsell                0.5.0    2018-06-12 [2] CRAN (R 4.0.0)
#  nlme                 * 3.1-147  2020-04-13 [3] CRAN (R 4.0.0)
#  nnet                   7.3-14   2020-04-26 [3] CRAN (R 4.0.0)
#  openssl                1.4.1    2019-07-18 [2] CRAN (R 4.0.0)
#  pillar                 1.4.4    2020-05-05 [2] CRAN (R 4.0.0)
#  pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.0.0)
#  plyr                   1.8.6    2020-03-03 [2] CRAN (R 4.0.0)
#  png                    0.1-7    2013-12-03 [2] CRAN (R 4.0.0)
#  prettyunits            1.1.1    2020-01-24 [2] CRAN (R 4.0.0)
#  progress               1.2.2    2019-05-16 [2] CRAN (R 4.0.0)
#  purrr                  0.3.4    2020-04-17 [1] CRAN (R 4.0.0)
#  qvalue                 2.20.0   2020-04-27 [2] Bioconductor
#  R6                     2.4.1    2019-11-12 [2] CRAN (R 4.0.0)
#  rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.0.0)
#  rappdirs               0.3.1    2016-03-28 [2] CRAN (R 4.0.0)
#  RColorBrewer         * 1.1-2    2014-12-07 [2] CRAN (R 4.0.0)
#  Rcpp                   1.0.4.6  2020-04-09 [1] CRAN (R 4.0.0)
#  RCurl                  1.98-1.2 2020-04-18 [2] CRAN (R 4.0.0)
#  readr                  1.3.1    2018-12-21 [2] CRAN (R 4.0.0)
#  recount                1.14.0   2020-04-27 [2] Bioconductor
#  rentrez                1.2.2    2019-05-02 [2] CRAN (R 4.0.0)
#  reshape2               1.4.4    2020-04-09 [2] CRAN (R 4.0.0)
#  rlang                  0.4.6    2020-05-02 [1] CRAN (R 4.0.0)
#  rngtools               1.5      2020-01-23 [2] CRAN (R 4.0.0)
#  rpart                  4.1-15   2019-04-12 [3] CRAN (R 4.0.0)
#  rprojroot              1.3-2    2018-01-03 [2] CRAN (R 4.0.0)
#  Rsamtools              2.4.0    2020-04-27 [2] Bioconductor
#  RSQLite                2.2.0    2020-01-07 [2] CRAN (R 4.0.0)
#  rstudioapi             0.11     2020-02-07 [2] CRAN (R 4.0.0)
#  rtracklayer            1.48.0   2020-04-27 [2] Bioconductor
#  S4Vectors            * 0.26.1   2020-05-16 [2] Bioconductor
#  scales                 1.1.1    2020-05-11 [2] CRAN (R 4.0.0)
#  segmented              1.1-0    2019-12-10 [2] CRAN (R 4.0.0)
#  sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.0.0)
#  stringi                1.4.6    2020-02-17 [2] CRAN (R 4.0.0)
#  stringr                1.4.0    2019-02-10 [2] CRAN (R 4.0.0)
#  SummarizedExperiment * 1.18.1   2020-04-30 [2] Bioconductor
#  survival               3.1-12   2020-04-10 [3] CRAN (R 4.0.0)
#  sva                  * 3.36.0   2020-04-27 [2] Bioconductor
#  tibble                 3.0.1    2020-04-20 [1] CRAN (R 4.0.0)
#  tidyr                  1.0.3    2020-05-07 [2] CRAN (R 4.0.0)
#  tidyselect             1.1.0    2020-05-11 [2] CRAN (R 4.0.0)
#  VariantAnnotation      1.34.0   2020-04-27 [2] Bioconductor
#  vctrs                  0.3.0    2020-05-11 [1] CRAN (R 4.0.0)
#  withr                  2.2.0    2020-04-20 [1] CRAN (R 4.0.0)
#  xfun                   0.13     2020-04-13 [2] CRAN (R 4.0.0)
#  XML                    3.99-0.3 2020-01-20 [2] CRAN (R 4.0.0)
#  xml2                   1.3.2    2020-04-23 [2] CRAN (R 4.0.0)
#  xtable                 1.8-4    2019-04-21 [2] CRAN (R 4.0.0)
#  XVector                0.28.0   2020-04-27 [2] Bioconductor
#  zlibbioc               1.34.0   2020-04-27 [2] Bioconductor
#
# [1] /users/lhuuki/R/4.0
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0/R/4.0/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0/R/4.0/lib64/R/library
#
#
