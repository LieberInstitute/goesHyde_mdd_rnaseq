library("data.table")
library("dplyr")
library("rio")
"%+%" <- function(a,b) paste(a, b, sep = "_")
"%&%" <- function(a,b) paste(a, b, sep = "")
data_dir <- "/dcl01/lieber/ImagingGenetics/QChen/Projects/RNAseq/DLPFC/HG38/CAUC_NC_SZ_MDD_BP/Data/"


##gene annotation
load(data_dir %&% "gene_annotation.DLPFC.CAUC.NC_SZ_MDD_BP.age13.rda")
genocode <- data.frame(pos_gene$Chr, pos_gene$name, pos_gene$Start, pos_gene$End, pos_gene$GeneType, pos_gene$Symbol)
setnames(genocode, old = c("pos_gene.Chr", "pos_gene.name", "pos_gene.Start", "pos_gene.End", "pos_gene.GeneType", "pos_gene.Symbol"),
         new = c("chr", "gene_id", "start", "end", "gene_type", "gene_name"))
export(genocode, "genocode.gtf", format = "gtf")


##snp annotation
load(data_dir %&% "snpMap.DLPFC.CAUC.NC_SZ_MDD_BP.age13.biallelic_snps.rda")
##care about only SNV type
snpMapNew <- snpMap[snpMap$Type == "SNV", ]
snpMapNew$VariantID <- substring(snpMapNew$chr_hg38, 4, 5) %+% snpMapNew$pos_hg38 %+% snpMapNew$newCount %+% snpMapNew$newRef %+% "hg38"
snpMapNew$snp_id_originalVCF <- "snp" %+% substring(snpMapNew$chr_hg38, 4, 5) %+% snpMapNew$pos_hg38
snpAnnot <- data.frame(substring(snpMapNew$chr_hg38, 4, 5), snpMapNew$pos_hg38, snpMapNew$VariantID, snpMapNew$newCount, snpMapNew$newRef, snpMapNew$snp_id_originalVCF, snpMapNew$name)
setnames(snpAnnot, old = c("substring.snpMapNew.chr_hg38..4..5.", "snpMapNew.pos_hg38", "snpMapNew.VariantID",
                           "snpMapNew.newCount", "snpMapNew.newRef", "snpMapNew.snp_id_originalVCF", "snpMapNew.name"),
         new = c("Chr", "Pos", "VariantID", "Ref_hg38", "Alt_hg38", "snp_id_originalVCF", "RSID"))
snpAnnot <- snpAnnot[snpAnnot$Pos != "2435598" & snpAnnot$Pos != "5916525", ]
head(snpAnnot[snpAnnot$Chr == 22, ])
write.table(snpAnnot, quote = FALSE, row.names = FALSE, file="snp_annotation.txt", sep = "\t")

##gene expression
load(data_dir %&% "log2_rpkm_gene.DLPFC.CAUC.NC_SZ_MDD_BP.age13.rda")
geneRpkm <- data.frame(log2_geneRpkm)
geneRpkm <- cbind(rownames(geneRpkm), data.frame(geneRpkm, row.names=NULL))
rownames(geneRpkm) <- c()
setnames(geneRpkm, old = "rownames(geneRpkm)", new = "TargetID")
write.table(geneRpkm, quote = FALSE, row.names = FALSE, file="gene_expression.txt", sep = "\t")

##covariates
load(data_dir %&% "covs_gene.DLPFC.CAUC.NC_SZ_MDD_BP.age13.rda")
write.table(covs_gene, quote = FALSE, row.names = FALSE, file="covariates.txt", sep="\t")


data_dir2 <- "/dcl01/lieber/ImagingGenetics/QChen/Projects/RNAseq/DLPFC/HG38/CAUC_NC_SZ_MDD_BP/Dosage/"

###mapping
load(data_dir %&% "pd.DLPFC.CAUC.NC_SZ_MDD_BP.age13.rda")
map <- data.frame(pd_output)


###genotype
for(i in 1:22) {
  load(data_dir2 %&% "postmortem_8runs.clinical_20170801_snps.chr" %&% i %&% ".rda")
  SNP$VariantID <- substring(SNP$chr_hg38, 4, 5) %+% SNP$pos_hg38 %+% SNP$A1 %+% SNP$A2 %+% "hg38"
  #SNP$VariantID <- SNP$pos_hg38
  gene <- data.frame(gen)
  snp <- data.frame(SNP)
  for(j in 1:392) {
    colnames(gene)[j] <- strsplit(colnames(gene)[j], "_")[[1]][2]
  }
  colnames(gene)[match(colnames(gene), map[, 1])] <- map[, 2][match(colnames(gene), map[, 1])]
  gene <- cbind(Id <- snp$VariantID, gene)
  setnames(gene, old = "Id <- snp$VariantID", new = "Id")
  write.table(gene, quote = FALSE, row.names = FALSE, file="genotype.chr" %&% i %&% ".txt", sep="\t")
}

dosage <- read.table("genotype.chr1.txt", header = TRUE)
for(i in 2:22) {
  dosage2 <- read.table("genotype.chr" %&% i %&% ".txt", header = TRUE)
  dosage <- rbind(dosage, dosage2)
}
dosage <- data.frame(dosage)
write.table(dosage, quote = FALSE, row.names = FALSE, file="genotype_dosage.txt", sep="\t")




###check match
genotype_snpID <- read.table("PredictDBPipeline/DLPFC_study/data/intermediate/genotypes/genotype.chr20.txt", header = T)
annotation_snpID <- read.table("PredictDBPipeline/DLPFC_study/data/intermediate/annotations/snp_annotation/snp_annot.chr20.txt", header = T)
before <- length(intersect(genotype_snpID$Id, annotation_snpID$varID))
annotation_snpID$newID <- annotation_snpID$chr %+% annotation_snpID$pos %+% annotation_snpID$effectAllele %+% annotation_snpID$refAllele %+% "hg38"
after <- length(intersect(genotype_snpID$Id, annotation_snpID$newID))

for(i in 1:nrow(genotype_snpID)) {
  genotype_snpID$pos[i] <- strsplit(as.character(genotype_snpID$Id[i]), "_")[[1]][2]
}

newMatch <- length(intersect(dosage$Id, annotation_snpID$varID))



expr <- readRDS("/users/ztong/PredictDBPipeline/DLPFC_study/data/intermediate/expression_phenotypes/gene_expression.RDS")
gene_annot <- readRDS("/users/ztong/PredictDBPipeline/DLPFC_study/data/intermediate/annotations/gene_annotation/gene_annot.RDS")
length(intersect(rownames(gene_annot), colnames(expr)))
dim(subset(gene_annot, gene_annot$chr == 2))

snp_annotation <- read.table("snp_annotation.txt", header = T)
gene_dosage <- read.table("genotype_dosage.txt", header = T)
length(intersect(gene_dosage$Id, snp_annotation$VariantID))
