####

### libraries
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)
library(sessioninfo)
library(here)

## load
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_exon.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_jxn.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_tx.Rdata"), verbose = TRUE)

## filter brain region
regInd <- which(colData(rse_gene)$BrainRegion == "sACC")
rse_gene <- rse_gene[, regInd]
rse_exon <- rse_exon[, regInd]
rse_jxn <- rse_jxn[, regInd]
rse_tx <- rse_tx[, regInd]

pd <- colData(rse_gene)

## load SNP data
load(here("eqtl", "genomewide", "rdas", "overlappingSNPs.rda"), verbose = TRUE) # snpMapKeep
load(here("genotype_data", "goesHyde_bipolarMdd_Genotypes_n593.rda"), verbose = TRUE)
head(snpMap)
corner(snp)
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


######################
# statistical model ##
######################
message(Sys.time(), " Get statistical model")
pd$PrimaryDx <- factor(pd$PrimaryDx,
    levels = c("Control", "Bipolar", "MDD")
)

mod <- model.matrix(~ PrimaryDx + Sex + as.matrix(mds[, 1:5]), data = pd)
colnames(mod)[grep("snpPC", colnames(mod))] <- colnames(mds)[1:5]


######################
# create SNP objects #
######################

theSnps <- SlicedData$new(as.matrix(snp))
theSnps$ResliceCombined(sliceSize = 50000)

snpspos <- snpMap[, c("SNP", "chr_hg38", "pos_hg38")]
colnames(snpspos) <- c("name", "chr", "pos")

# remove snpMap snp mds
rm(snpMap, snp, mds)

#####################
## calculate rpkm ##
####################
message(Sys.time(), " Get rpkm values")
geneRpkm <- recount::getRPKM(rse_gene, "Length")
exonRpkm <- recount::getRPKM(rse_exon, "Length")
rowData(rse_jxn)$Length <- 100
jxnRp10m <- recount::getRPKM(rse_jxn, "Length")
txTpm <- assays(rse_tx)$tpm


#######################
####### do PCA ########
#######################
message(Sys.time(), " Do PCA")

pca_rda_file <- here("eqtl", "genomewide", "rdas", "pcs_4features_sacc.rda")

if(!file.exists(pca_rda_file)){
    pcaGene <- prcomp(t(log2(geneRpkm + 1)))
    kGene <- num.sv(log2(geneRpkm + 1), mod)
    genePCs <- pcaGene$x[, 1:kGene]

    pcaExon <- prcomp(t(log2(exonRpkm + 1)))
    kExon <- num.sv(log2(exonRpkm + 1), mod, vfilter = 50000)
    exonPCs <- pcaExon$x[, 1:kExon]

    pcaJxn <- prcomp(t(log2(jxnRp10m + 1)))
    kJxn <- num.sv(log2(jxnRp10m + 1), mod, vfilter = 50000)
    jxnPCs <- pcaJxn$x[, 1:kJxn]

    pcaTx <- prcomp(t(log2(txTpm + 1)))
    kTx <- num.sv(log2(txTpm + 1), mod, vfilter = 50000)
    txPCs <- pcaTx$x[, 1:kTx]

    save(genePCs, exonPCs, jxnPCs, txPCs, file = pca_rda_file)

}else{
    load(pca_rda_file, verbose = TRUE)
}

covsGene <- SlicedData$new(t(cbind(mod[, -1], genePCs)))
covsExon <- SlicedData$new(t(cbind(mod[, -1], exonPCs)))
covsJxn <- SlicedData$new(t(cbind(mod[, -1], jxnPCs)))
covsTx <- SlicedData$new(t(cbind(mod[, -1], txPCs)))


##########################
### feature annotation ###
##########################

###### gene level
posGene <- as.data.frame(rowRanges(rse_gene))[, 1:3]
posGene$name <- rownames(posGene)
posGene <- posGene[, c(4, 1:3)]

##### exon level
posExon <- as.data.frame(rowRanges(rse_exon))[, 1:3]
posExon$name <- rownames(posExon)
posExon <- posExon[, c(4, 1:3)]

##### junction level
posJxn <- as.data.frame(rowRanges(rse_jxn))[, 1:3]
posJxn$name <- rownames(posJxn)
posJxn <- posJxn[, c(4, 1:3)]
names(posJxn)[2:4] <- c("Chr", "Start", "End")

##### transcript level
posTx <- as.data.frame(rowRanges(rse_tx))[, 1:3]
posTx$name <- rownames(posTx)
posTx <- posTx[, c(4, 1:3)]
names(posTx)[2:4] <- c("Chr", "Start", "End")

#############################
### sliced expression data ##
geneSlice <- SlicedData$new(log2(geneRpkm + 1))
exonSlice <- SlicedData$new(log2(exonRpkm + 1))
jxnSlice <- SlicedData$new(log2(jxnRp10m + 1))
txSlice <- SlicedData$new(log2(txTpm + 1))

geneSlice$ResliceCombined(sliceSize = 5000)
exonSlice$ResliceCombined(sliceSize = 5000)
jxnSlice$ResliceCombined(sliceSize = 5000)
txSlice$ResliceCombined(sliceSize = 5000)


##########################
### Run EQTLs ############
##########################
message(Sys.time(), " EQTLs")

print("Starting eQTLs")
# takes a long time
meGene_rda <- "matrixEqtl_output_sacc_genomewide_gene.rda"
if(!file.exists(meGene_rda)){
    meGene <- Matrix_eQTL_main(
        snps = theSnps, gene = geneSlice,
        cvrt = covsGene, output_file_name.cis = ".ctxt",
        pvOutputThreshold.cis = .1, pvOutputThreshold = 0,
        snpspos = snpspos, genepos = posGene,
        useModel = modelLINEAR, cisDist = 5e5,
        pvalue.hist = 100, min.pv.by.genesnp = TRUE
    )
    save(meGene, file = meGene_rda)
}else{
    load(meGene_rda, verbose = TRUE)
}

meExon_rda <- "matrixEqtl_output_sacc_genomewide_exon.rda"
if(!file.exists(meExon_rda)){
    meExon <- Matrix_eQTL_main(
        snps = theSnps, gene = exonSlice,
        cvrt = covsExon, output_file_name.cis = ".ctxt",
        pvOutputThreshold.cis = .1, pvOutputThreshold = 0,
        snpspos = snpspos, genepos = posExon,
        useModel = modelLINEAR, cisDist = 5e5,
        pvalue.hist = 100, min.pv.by.genesnp = TRUE
    )
    save(meExon, file = meExon_rda)
}else{
    load(meExon_rda, verbose = TRUE)
}

meJxn_rda <- "matrixEqtl_output_sacc_genomewide_jxn.rda"
if(!file.exists(meJxn_rda)){
    meJxn <- Matrix_eQTL_main(
        snps = theSnps, gene = jxnSlice,
        cvrt = covsJxn, output_file_name.cis = ".ctxt",
        pvOutputThreshold.cis = .1, pvOutputThreshold = 0,
        snpspos = snpspos, genepos = posJxn,
        useModel = modelLINEAR, cisDist = 5e5,
        pvalue.hist = 100, min.pv.by.genesnp = TRUE
    )
    save(meJxn, file = meJxn_rda)
}else{
    load(meJxn_rda, verbose = TRUE)
}

meTx_rda <- "matrixEqtl_output_sacc_genomewide_tx.rda"
if(!file.exists(meTx_rda)){
    meTx <- Matrix_eQTL_main(
        snps = theSnps, gene = txSlice,
        cvrt = covsTx, output_file_name.cis = ".ctxt",
        pvOutputThreshold.cis = .1, pvOutputThreshold = 0,
        snpspos = snpspos, genepos = posTx,
        useModel = modelLINEAR, cisDist = 5e5,
        pvalue.hist = 100, min.pv.by.genesnp = TRUE
    )
    save(meTx, file = meTx_rda)
}else{
    load(meTx_rda, verbose = TRUE)
}


######################
###### annotate ######
#####################

message(Sys.time(), " Annotate and Save")
# extract
geneEqtl <- meGene$cis$eqtls
exonEqtl <- meExon$cis$eqtls
jxnEqtl <- meJxn$cis$eqtls
txEqtl <- meTx$cis$eqtls

################################
# add gene annotation info #####
################################


geneEqtl$Symbol <- rowRanges(rse_gene)$Symbol[match(geneEqtl$gene, rownames(rse_gene))]
geneEqtl$EnsemblGeneID <- rowRanges(rse_gene)$ensemblID[match(geneEqtl$gene, rownames(rse_gene))]
geneEqtl$Type <- "Gene"
geneEqtl$Class <- "InGen" # in gencode
geneEqtl <- DataFrame(geneEqtl)
# geneEqtl$gene_type = rowRanges(rse_gene)$gene_type[match(geneEqtl$gene, rownames(rse_gene))]

exonEqtl$Symbol <- rowRanges(rse_exon)$Symbol[match(exonEqtl$gene, rownames(rse_exon))]
exonEqtl$EnsemblGeneID <- rowRanges(rse_exon)$ensemblID[match(exonEqtl$gene, rownames(rse_exon))]
exonEqtl$Type <- "Exon"
exonEqtl$Class <- "InGen"
exonEqtl <- DataFrame(exonEqtl)
# exonEqtl$gene_type = rowRanges(rse_exon)$gene_type[match(exonEqtl$gene, rownames(rse_exon))]

jxnEqtl$Symbol <- rowRanges(rse_jxn)$newGeneSymbol[match(jxnEqtl$gene, rownames(rse_jxn))]
jxnEqtl$EnsemblGeneID <- rowRanges(rse_jxn)$newGeneID[match(jxnEqtl$gene, rownames(rse_jxn))]
jxnEqtl$Type <- "Jxn"
jxnEqtl$Class <- rowRanges(rse_jxn)$Class[match(jxnEqtl$gene, rownames(rse_jxn))]
jxnEqtl <- DataFrame(jxnEqtl)
# jxnEqtl$gene_type = rowRanges(rse_jxn)$gene_type[match(jxnEqtl$gene, rownames(rse_jxn))]

txEqtl$Symbol <- rowRanges(rse_tx)$gene_name[match(txEqtl$gene, rownames(rse_tx))]
txEqtl$EnsemblGeneID <- ss(rowRanges(rse_tx)$gene_id[match(txEqtl$gene, rownames(rse_tx))], "\\.", 1)
txEqtl$Type <- "Tx"
txEqtl$Class <- "InGen"
txEqtl <- DataFrame(txEqtl)
# txEqtl$gene_type = rowRanges(rse_tx)$gene_type[match(txEqtl$gene, rownames(rse_tx))]

# remove rse objects
rm(rse_gene, rse_exon, rse_jxn, rse_tx)

# merge
allEqtl <- rbind(geneEqtl, exonEqtl, jxnEqtl, txEqtl)
rm(geneEqtl, exonEqtl, jxnEqtl, txEqtl)
# allEqtl$gencodeTx <- CharacterList(c(
#     as.list(rowRanges(rse_gene)$gencodeTx[match(
#         geneEqtl$gene,
#         rownames(rse_gene)
#     )]),
#     as.list(rowRanges(rse_exon)$gencodeTx[match(exonEqtl$gene, rownames(rse_exon))]),
#     as.list(rowRanges(rse_jxn)$gencodeTx[match(jxnEqtl$gene, rownames(rse_jxn))]),
#     as.list(txEqtl$gene)
# ))
save(allEqtl, file = "mergedEqtl_output_sacc_genomewide_4features.rda", compress = TRUE)


allEqtlFDR01 <- allEqtl[which(allEqtl$FDR < 0.01), ]
save(allEqtlFDR01, file = "mergedEqtl_output_sacc_genomewide_4features_FDR01.rda", compress = TRUE)

# sgejobs::job_single("genomewide_run_eqtls_sacc", memory = "150G",create_shell = TRUE, command = "Rscript genomewide_run_eqtls_sacc.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# [1] "2020-06-01 14:42:39 EDT"
#      user    system   elapsed
# 12311.991   464.563 12819.941
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
#  date     2020-06-01
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
#  RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.0.0)
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
# **** Job ends ****
# Mon Jun  1 14:42:45 EDT 2020
#
