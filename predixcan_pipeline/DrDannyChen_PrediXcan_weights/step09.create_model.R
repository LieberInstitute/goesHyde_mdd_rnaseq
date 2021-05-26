argv <- commandArgs(trailingOnly = TRUE)
source("GTEx_Tissue_Wide_CV_elasticNet.R")

tis <- argv[1]
chrom <- argv[2]
alpha <- as.numeric(argv[3])
window <- as.numeric(argv[4])

data_dir <- "../data/intermediate/"

expression_RDS <- data_dir %&% "expression_phenotypes/" %&% "gene_expression.RDS"
geno_file <- data_dir %&% "genotypes/" %&% "genotype.chr" %&% chrom %&% ".txt"
gene_annot_RDS <- data_dir %&% "annotations/gene_annotation/gene_annot.RDS"
snp_annot_RDS <- data_dir %&% "annotations/snp_annotation/snp_annot.chr" %&% chrom %&% ".RDS"
n_k_folds <- 10
out_dir <- data_dir %&% "model_by_chr/"
snpset <- "1KG"

TW_CV_model(expression_RDS, geno_file, gene_annot_RDS, snp_annot_RDS,
    n_k_folds, alpha, out_dir, tis, chrom, snpset, window)
