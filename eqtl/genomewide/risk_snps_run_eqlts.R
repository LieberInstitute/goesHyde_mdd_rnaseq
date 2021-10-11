
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

## Risk SNPS 
## HG19
risk_snps <- read.csv(here("eqtl","genomewide","pgc3-supptable2.csv"))
risk_snps$BP <- as.numeric(gsub(",","",risk_snps$BP))
risk_snps$pos_hg19 <- paste0(risk_snps$CHR, ":", risk_snps$BP)
nrow(risk_snps)
#64

## Match ?
load(here("genotype_data","dbSNP.Rdata"), verbose = TRUE)

### match to SNPs
risk_snpsGR = GRanges(paste0("chr", risk_snps$CHR), 
                      IRanges(start = risk_snps$BP,width=1,names=risk_snps$SNP))
oo = findOverlaps(risk_snpsGR, dbSnp142)

risk_snps$rsNumGuess[queryHits(oo)] = dbSnp142$RefSNP_id[subjectHits(oo)]

head(risk_snps)
table(risk_snps$SNP == risk_snps$rsNumGuess)
rm(dbSnp142, dbSnp149)

## load SNP data
load(here("genotype_data", "goesHyde_bipolarMdd_Genotypes.rda"), verbose = TRUE)
rownames(snp) <- rownames(snpMap) <- snpMap$SNP
snpMap$pos_hg19 <- paste0(snpMap$CHR, ":", snpMap$Pos)


intersect(risk_snps$SNP, snpMap$rsNumGuess)
intersect(risk_snps$pos_hg19, snpMap$pos_hg19)

