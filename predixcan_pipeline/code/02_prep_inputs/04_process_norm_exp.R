library(GenomicRanges)

norm_exp <-
  read.delim(
    text = gsub(
      "#",
      "",
      readLines(
        "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/goesHyde_mdd_rnaseq_Amygdala.expression.bed.gz"
      )
    ),
    header = TRUE,
    check.names = FALSE,
    sep = ""
  )

norm_exp = norm_exp[-c(1:3)]

colnames(norm_exp)[1] <- "Gene_Name"

# Transpose the data to enable you rub peer factors
n = norm_exp$Gene_Name
norm_exp_transpose <- as.data.frame(t(norm_exp[, -1]))
colnames(norm_exp_transpose) <- n

write.table(
  norm_exp_transpose,
  file = "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/02_prep_inputs/transformed_expression.txt",
  sep = "\t",
  row.names = TRUE
)

# peer_factors = read.csv(file = "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/goesHyde_mdd_rnaseq_Amygdala.combined_covariates.txt", header = TRUE, sep = "\t", row.names = 1)
# colnames(peer_factors) = rownames(norm_exp_transpose)

# write out a covariates matrix
# write.table(peer_factors, file = "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/02_prep_inputs/covariates_matrix.txt", sep = "\t",
#             row.names = TRUE)

# expression = gene_exp_transpose
