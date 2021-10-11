
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)
library(sessioninfo)
library(here)
library(VariantAnnotation)

## load
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)

## filter brain region
rse_gene <- rse_gene[, rse_gene$BrainRegion == "Amygdala"]
dim(rse_gene)
pd <- colData(rse_gene)

## load SNP data
load(here("genotype_data", "goesHyde_bipolarMdd_Genotypes.rda"), verbose = TRUE) # mds snp snpMap
head(snpMap)
corner(snp)
dim(snp)
# [1] 8892060     595


snpMap$pos_hg19 <- paste0(snpMap$CHR, ":", snpMap$POS)
rownames(snp) <- rownames(snpMap) <- snpMap$SNP

# table(rownames(snpMap) %in% rownames(snpMapKeep))
# snpInd <- which(rownames(snpMap) %in% rownames(snpMapKeep) & !is.na(snpMap$pos_hg38))

## filter snps w/ no HG38 pos
has_pos <- !is.na(snpMap$pos_hg38)
table(has_pos)
# FALSE    TRUE 
# 92045 8800015 

snpMap <- snpMap[has_pos, ]
snp <- snp[has_pos, ]

dim(snp)
# [1] 8401829     540

#####################
# filter brain region
# make mds and snp dimensions equal to N
# (repeat rows or columns for BrNum replicates)

all(pd$BrNum %in% colnames(snp))
all(pd$BrNum %in% rownames(mds))
all(pd$BrNum %in% rownames(mds))

mds <- mds[pd$BrNum, ]
snp <- snp[, pd$BrNum]
# rownames(mds) <- colnames(snp) <- pd$RNum
# snpMap$maf <- rowSums(snp, na.rm = TRUE) / (2 * rowSums(!is.na(snp)))
# summary(snpMap$maf)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.005556 0.037963 0.124074 0.171798 0.288889 0.515741

######################
# statistical model ##
######################
message(Sys.time(), " Get statistical model")
pd$PrimaryDx <- factor(pd$PrimaryDx,
    levels = c("Control", "Bipolar", "MDD")
)

mod <- cbind(model.matrix(~ PrimaryDx + Sex, data = pd), as.matrix(mds[, 1:5]))
colnames(mod)

######################
# create SNP objects #
######################
message(Sys.time(), " Create SlicedData")
theSnps <- SlicedData$new(as.matrix(snp))
theSnps$ResliceCombined(sliceSize = 50000)

snpspos <- snpMap[, c("SNP", "chr_hg38", "pos_hg38")]
colnames(snpspos) <- c("name", "chr", "pos")

# remove snpMap snp mds
rm(snpMap, snp, mds)

#### calculate rpkm ####
message(Sys.time(), " Get rpkm values")
geneRpkm <- recount::getRPKM(rse_gene, "Length")

#### Load PCA ####
message(Sys.time(), " Load PCA")

load(here("eqtl", "data", "featuresPCs","genePCsAmygdala.rda"), verbose = TRUE)
covsGene <- SlicedData$new(t(cbind(mod[, -1], PCs)))


#### feature annotation ####
## gene level
posGene <- as.data.frame(rowRanges(rse_gene))[, 1:3]
posGene$name <- rownames(posGene)
posGene <- posGene[, c(4, 1:3)]


#### sliced expression data ###
geneSlice <- SlicedData$new(log2(geneRpkm + 1))
geneSlice$ResliceCombined(sliceSize = 5000)

### Run EQTLs ####
message(Sys.time(), " EQTLs")

# takes a long time
meGene_rda <- here("eqtl", "data", "matrixEQTL_out", "matrixEqtl_output_amyg_genomewide_gene.rda")

meGene <- Matrix_eQTL_main(
    snps = theSnps, gene = geneSlice,
    cvrt = covsGene, output_file_name.cis = ".ctxt",
    pvOutputThreshold.cis = .1, pvOutputThreshold = 0,
    snpspos = snpspos, genepos = posGene,
    useModel = modelLINEAR, cisDist = 5e5,
    pvalue.hist = 100, min.pv.by.genesnp = TRUE
)
save(meGene, file = meGene_rda)

#### annotate ####
message(Sys.time(), " Annotate and Save")
# extract
geneEqtl <- meGene$cis$eqtls

## add gene annotation info 
geneEqtl$Symbol <- rowRanges(rse_gene)$Symbol[match(geneEqtl$gene, rownames(rse_gene))]
geneEqtl$EnsemblGeneID <- rowRanges(rse_gene)$ensemblID[match(geneEqtl$gene, rownames(rse_gene))]
geneEqtl$Type <- "Gene"
geneEqtl$Class <- "InGen" # in gencode
geneEqtl <- DataFrame(geneEqtl)
# geneEqtl$gene_type = rowRanges(rse_gene)$gene_type[match(geneEqtl$gene, rownames(rse_gene))]

save(geneEqtl, file = here("eqtl", "data", "matrixEQTL_out", "matrixEqtl_output_amyg_genomewide_gene_annotate.rda")
)
# sgejobs::job_single("matrixEQTL_genomewide", memory = "150G",create_shell = TRUE, command = "Rscript matrixEQTL_genomewide.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
