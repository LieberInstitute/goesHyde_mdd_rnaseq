## Partially adapted from clean_data.R
library('SummarizedExperiment')
library('sessioninfo')

## Load data
load('degradation_rse_MDDseq_BothRegions.rda', verbose = TRUE)
rse_mdd <- cov_rse
rse_mdd$Experiment <- "GoesMDD"
load('/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/data/degradation_rse_BipSeq_BothRegions.rda', verbose = TRUE)
rse_bip <- cov_rse
rse_bip$Experiment <- "ZandiBPD"

dim(rse_mdd)
# [1] 1000  632
dim(rse_bip)
# [1] 1000  511

## Are the degradation regions the same ones? Yes
stopifnot(identical(rowRanges(rse_mdd), rowRanges(rse_bip)))

## Do the samples overlap? 9 do as Emily had noted already in clean_data.R
table(rse_mdd$SAMPLE_ID %in% rse_bip$RNum)
# FALSE  TRUE
#   623     9
table(rse_bip$RNum %in% rse_mdd$SAMPLE_ID)
# FALSE  TRUE
#   502     9

## 9 Control samples also got used in Goes -> drop from Zandi so there's no duplicates
rse_bip <- rse_bip[, -which(rse_bip$RNum %in% colnames(rse_mdd))]	# 502

## combine
# make colData consistent
rse_mdd$AgeDeath <- rse_mdd$Age
rse_mdd$RNum <- rse_mdd$SAMPLE_ID
rse_mdd$BrNum <- as.character(rse_mdd$BrNum)
colKeep <- c("RNum","BrNum","Sex","Race","AgeDeath","BrainRegion","RIN","PrimaryDx","overallMapRate","totalAssignedGene","mitoRate","rRNA_rate","Experiment")
colData(rse_mdd) <- colData(rse_mdd)[,colKeep]
colData(rse_bip) <- colData(rse_bip)[,colKeep]

rse_both <- cbind(rse_mdd, rse_bip)

## drop
qc <- read.csv("../qc_checks/qc_dropping_results.csv", stringsAsFactors = FALSE)
qc <- qc[rowSums(qc[,13:16])>0,]
qc <- rbind(qc, "R17779")	# Tiny # of reads
rse_both <- rse_both[,-which(rse_both$RNum %in% qc$SAMPLE_ID | rse_both$RNum %in% c("R17538","R18853") |
						rse_both$PrimaryDx == "Other" |
						rse_both$overallMapRate <0.5) ]

rse_both$PrimaryDx <- droplevels(rse_both$PrimaryDx)
rse_both$PrimaryDx <- relevel(rse_both$PrimaryDx, ref="MDD")


## Compare against the merged one from before
load('rse_gene_GoesZandi_n1093.rda', verbose = TRUE)

## They are not identical...
identical(colData(rse_gene), colData(rse_both))

stopifnot(identical(colnames(colData(rse_gene)), colnames(colData(rse_both))))
stopifnot(all(sapply(seq_along(colData(rse_gene)), function(i) {
    identical(colData(rse_gene)[, i], colData(rse_both)[, i])
})))
stopifnot(identical(colnames(rse_gene), colnames(rse_both)))

## But they are equivalent, hm...
## Oh well, now we have a degradation object we can use that matches
## the rse_gene MDD + BiP seq merged one Emily made
testthat::expect_equivalent(colData(rse_gene), colData(rse_both))

## Rename to the convention we've used
cov_rse <- rse_both
save(cov_rse, file = 'degradation_rse_MDDseq_BiPSeq_BothRegions.Rdata')


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# sgejobs::job_single('merge_MDDseq_BiPseq_degradation', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript merge_MDDseq_BiPseq_degradation.R")

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.6.1 Patched (2019-10-31 r77350)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2020-04-02
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version  date       lib source
#  acepack                1.4.1    2016-10-29 [2] CRAN (R 3.6.1)
#  AnnotationDbi          1.48.0   2019-10-29 [1] Bioconductor
#  askpass                1.1      2019-01-13 [1] CRAN (R 3.6.1)
#  assertthat             0.2.1    2019-03-21 [2] CRAN (R 3.6.1)
#  backports              1.1.5    2019-10-02 [1] CRAN (R 3.6.1)
#  base64enc              0.1-3    2015-07-28 [2] CRAN (R 3.6.1)
#  Biobase              * 2.46.0   2019-10-29 [2] Bioconductor
#  BiocFileCache          1.10.2   2019-11-08 [1] Bioconductor
#  BiocGenerics         * 0.32.0   2019-10-29 [1] Bioconductor
#  BiocParallel         * 1.20.1   2019-12-21 [1] Bioconductor
#  biomaRt                2.42.0   2019-10-29 [1] Bioconductor
#  Biostrings             2.54.0   2019-10-29 [1] Bioconductor
#  bit                    1.1-15.2 2020-02-10 [2] CRAN (R 3.6.1)
#  bit64                  0.9-7    2017-05-08 [2] CRAN (R 3.6.1)
#  bitops                 1.0-6    2013-08-17 [2] CRAN (R 3.6.1)
#  blob                   1.2.1    2020-01-20 [2] CRAN (R 3.6.1)
#  BSgenome               1.54.0   2019-10-29 [1] Bioconductor
#  bumphunter             1.28.0   2019-10-29 [1] Bioconductor
#  checkmate              2.0.0    2020-02-06 [1] CRAN (R 3.6.1)
#  cli                    2.0.2    2020-02-28 [1] CRAN (R 3.6.1)
#  cluster                2.1.0    2019-06-19 [3] CRAN (R 3.6.1)
#  codetools              0.2-16   2018-12-24 [3] CRAN (R 3.6.1)
#  colorout             * 1.2-2    2019-10-31 [1] Github (jalvesaq/colorout@641ed38)
#  colorspace             1.4-1    2019-03-18 [2] CRAN (R 3.6.1)
#  crayon                 1.3.4    2017-09-16 [1] CRAN (R 3.6.1)
#  curl                   4.3      2019-12-02 [1] CRAN (R 3.6.1)
#  data.table             1.12.8   2019-12-09 [1] CRAN (R 3.6.1)
#  DBI                    1.1.0    2019-12-15 [2] CRAN (R 3.6.1)
#  dbplyr                 1.4.2    2019-06-17 [1] CRAN (R 3.6.1)
#  DelayedArray         * 0.12.2   2020-01-06 [2] Bioconductor
#  derfinder              1.20.0   2019-10-29 [1] Bioconductor
#  derfinderHelper        1.20.0   2019-10-29 [1] Bioconductor
#  digest                 0.6.25   2020-02-23 [1] CRAN (R 3.6.1)
#  doRNG                  1.8.2    2020-01-27 [2] CRAN (R 3.6.1)
#  downloader             0.4      2015-07-09 [1] CRAN (R 3.6.1)
#  dplyr                  0.8.5    2020-03-07 [1] CRAN (R 3.6.1)
#  fansi                  0.4.1    2020-01-08 [1] CRAN (R 3.6.1)
#  foreach                1.4.8    2020-02-09 [2] CRAN (R 3.6.1)
#  foreign                0.8-72   2019-08-02 [3] CRAN (R 3.6.1)
#  Formula                1.2-3    2018-05-03 [2] CRAN (R 3.6.1)
#  GenomeInfoDb         * 1.22.0   2019-10-29 [1] Bioconductor
#  GenomeInfoDbData       1.2.2    2019-10-28 [2] Bioconductor
#  GenomicAlignments      1.22.1   2019-11-12 [1] Bioconductor
#  GenomicFeatures        1.38.2   2020-02-15 [1] Bioconductor
#  GenomicFiles           1.22.0   2019-10-29 [1] Bioconductor
#  GenomicRanges        * 1.38.0   2019-10-29 [1] Bioconductor
#  GEOquery               2.54.1   2019-11-18 [1] Bioconductor
#  ggplot2                3.3.0    2020-03-05 [1] CRAN (R 3.6.1)
#  glue                   1.3.2    2020-03-12 [1] CRAN (R 3.6.1)
#  gridExtra              2.3      2017-09-09 [2] CRAN (R 3.6.1)
#  gtable                 0.3.0    2019-03-25 [2] CRAN (R 3.6.1)
#  Hmisc                  4.3-1    2020-02-07 [1] CRAN (R 3.6.1)
#  hms                    0.5.3    2020-01-08 [2] CRAN (R 3.6.1)
#  htmlTable              1.13.3   2019-12-04 [1] CRAN (R 3.6.1)
#  htmltools              0.4.0    2019-10-04 [1] CRAN (R 3.6.1)
#  htmlwidgets            1.5.1    2019-10-08 [1] CRAN (R 3.6.1)
#  httpuv                 1.5.2    2019-09-11 [1] CRAN (R 3.6.1)
#  httr                   1.4.1    2019-08-05 [1] CRAN (R 3.6.1)
#  IRanges              * 2.20.2   2020-01-13 [1] Bioconductor
#  iterators              1.0.12   2019-07-26 [2] CRAN (R 3.6.1)
#  jpeg                   0.1-8.1  2019-10-24 [2] CRAN (R 3.6.1)
#  jsonlite               1.6.1    2020-02-02 [2] CRAN (R 3.6.1)
#  knitr                  1.28     2020-02-06 [1] CRAN (R 3.6.1)
#  later                  1.0.0    2019-10-04 [1] CRAN (R 3.6.1)
#  lattice                0.20-38  2018-11-04 [3] CRAN (R 3.6.1)
#  latticeExtra           0.6-29   2019-12-19 [2] CRAN (R 3.6.1)
#  lifecycle              0.2.0    2020-03-06 [1] CRAN (R 3.6.1)
#  limma                  3.42.2   2020-02-03 [1] Bioconductor
#  locfit                 1.5-9.1  2013-04-20 [2] CRAN (R 3.6.1)
#  magrittr               1.5      2014-11-22 [1] CRAN (R 3.6.1)
#  Matrix                 1.2-17   2019-03-22 [3] CRAN (R 3.6.1)
#  matrixStats          * 0.56.0   2020-03-13 [1] CRAN (R 3.6.1)
#  memoise                1.1.0    2017-04-21 [2] CRAN (R 3.6.1)
#  munsell                0.5.0    2018-06-12 [2] CRAN (R 3.6.1)
#  nnet                   7.3-12   2016-02-02 [3] CRAN (R 3.6.1)
#  openssl                1.4.1    2019-07-18 [1] CRAN (R 3.6.1)
#  pillar                 1.4.3    2019-12-20 [1] CRAN (R 3.6.1)
#  pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 3.6.1)
#  plyr                   1.8.5    2019-12-10 [2] CRAN (R 3.6.1)
#  png                    0.1-7    2013-12-03 [2] CRAN (R 3.6.1)
#  prettyunits            1.1.1    2020-01-24 [1] CRAN (R 3.6.1)
#  progress               1.2.2    2019-05-16 [1] CRAN (R 3.6.1)
#  promises               1.1.0    2019-10-04 [1] CRAN (R 3.6.1)
#  purrr                  0.3.3    2019-10-18 [2] CRAN (R 3.6.1)
#  qvalue                 2.18.0   2019-10-29 [1] Bioconductor
#  R6                     2.4.1    2019-11-12 [2] CRAN (R 3.6.1)
#  rappdirs               0.3.1    2016-03-28 [1] CRAN (R 3.6.1)
#  RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 3.6.1)
#  Rcpp                   1.0.3    2019-11-08 [1] CRAN (R 3.6.1)
#  RCurl                  1.98-1.1 2020-01-19 [2] CRAN (R 3.6.1)
#  readr                  1.3.1    2018-12-21 [1] CRAN (R 3.6.1)
#  recount                1.12.1   2019-11-06 [1] Bioconductor
#  rentrez                1.2.2    2019-05-02 [1] CRAN (R 3.6.1)
#  reshape2               1.4.3    2017-12-11 [2] CRAN (R 3.6.1)
#  rlang                  0.4.5    2020-03-01 [1] CRAN (R 3.6.1)
#  rmote                * 0.3.4    2019-10-31 [1] Github (cloudyr/rmote@fbce611)
#  rngtools               1.5      2020-01-23 [2] CRAN (R 3.6.1)
#  rpart                  4.1-15   2019-04-12 [3] CRAN (R 3.6.1)
#  Rsamtools              2.2.3    2020-02-23 [1] Bioconductor
#  RSQLite                2.2.0    2020-01-07 [2] CRAN (R 3.6.1)
#  rstudioapi             0.11     2020-02-07 [2] CRAN (R 3.6.1)
#  rtracklayer            1.46.0   2019-10-29 [1] Bioconductor
#  S4Vectors            * 0.24.3   2020-01-18 [1] Bioconductor
#  scales                 1.1.0    2019-11-18 [2] CRAN (R 3.6.1)
#  servr                  0.16     2020-03-02 [1] CRAN (R 3.6.1)
#  sessioninfo          * 1.1.1    2018-11-05 [1] CRAN (R 3.6.1)
#  stringi                1.4.6    2020-02-17 [2] CRAN (R 3.6.1)
#  stringr                1.4.0    2019-02-10 [1] CRAN (R 3.6.1)
#  SummarizedExperiment * 1.16.1   2019-12-19 [1] Bioconductor
#  survival               3.1-11   2020-03-07 [1] CRAN (R 3.6.1)
#  testthat               2.3.2    2020-03-02 [1] CRAN (R 3.6.1)
#  tibble                 2.1.3    2019-06-06 [1] CRAN (R 3.6.1)
#  tidyr                  1.0.2    2020-01-24 [1] CRAN (R 3.6.1)
#  tidyselect             1.0.0    2020-01-27 [2] CRAN (R 3.6.1)
#  VariantAnnotation      1.32.0   2019-10-29 [1] Bioconductor
#  vctrs                  0.2.4    2020-03-10 [1] CRAN (R 3.6.1)
#  withr                  2.1.2    2018-03-15 [2] CRAN (R 3.6.1)
#  xfun                   0.12     2020-01-13 [1] CRAN (R 3.6.1)
#  XML                    3.99-0.3 2020-01-20 [2] CRAN (R 3.6.1)
#  xml2                   1.2.2    2019-08-09 [2] CRAN (R 3.6.1)
#  XVector                0.26.0   2019-10-29 [1] Bioconductor
#  zlibbioc               1.32.0   2019-10-29 [2] Bioconductor
#
# [1] /users/lcollado/R/3.6.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/library
