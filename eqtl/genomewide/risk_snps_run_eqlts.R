
library(SummarizedExperiment)
library(jaffelab)
# library(MatrixEQTL)
library(sva)
library(sessioninfo)
library(here)
library(dplyr)

## load ## loop eventually ?
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
# load(here("exprs_cutoff", "rse_exon.Rdata"), verbose = TRUE)
# load(here("exprs_cutoff", "rse_jxn.Rdata"), verbose = TRUE)
# load(here("exprs_cutoff", "rse_tx.Rdata"), verbose = TRUE)

## load SNP data
load(here("genotype_data", "goesHyde_bipolarMdd_Genotypes.rda"), verbose = TRUE)
rownames(snp) <- rownames(snpMap) <- snpMap$SNP
snpMap$pos_hg19 <- paste0(snpMap$CHR, ":", snpMap$POS)

risk_snps <- read.csv(here("eqtl","genomewide","pgc3-supptable2.csv"))
risk_snps$BP <- as.numeric(gsub(",","",risk_snps$BP))
risk_snps$pos_hg19 <- paste0(risk_snps$CHR, ":", risk_snps$BP)
nrow(risk_snps)
#64

intersect(risk_snps$pos_hg19, snpMap$pos_hg19)

head(risk_snps$BP)
head(snpMap$POS)


## filter snps
has_pos <- !is.na(snpMap$pos_hg38)
snpMap <- snpMap[has_pos, ]
snp <- snp[has_pos, ]



