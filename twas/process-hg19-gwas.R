library(data.table)
library(dplyr)

setDTthreads(1)

# Load snpMap
load(
    "/dcl01/lieber/ajaffe/Brain/Imputation/Subj_Cleaned/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38.Rdata"
)

snpMap <- as.data.table(snpMap)

# Read in Fernando's GWAS
hg19_gwas <-
    fread(
        "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/PGC_MDD_sumstats/PGC_UKB_23andMe_depression_genome-wide.txt"
    )

table(snpMap$rsNumGuess %in% hg19_gwas$MarkerName)
table(snpMap$name %in% hg19_gwas$MarkerName)
# > nrow(hg19_gwas)
# [1] 8098588
# > 4408597/ nrow(hg19_gwas)
# [1] 0.5443661
# > 4408597+6579814
# [1] 10988411
# > 4408597 / nrow(snpMap)
# [1] 0.4012042
# > 4410625 / nrow(snpMap)
# [1] 0.4013888

snpMap <- snpMap[snpMap$name %in% hg19_gwas$MarkerName, ]

# Merge in the hg38 coordinates based on hg19 coordinates
hg38_gwas <-
    merge(hg19_gwas, snpMap[, c("chr_hg38", "pos_hg38", "name")], by.x = "MarkerName", by.y = "name")

# Drop hg19 coords from gwas
# hg38_gwas <- hg38_gwas[, -c("hg19_key", "CHR", "BP")]

# reorder columns
col_order <- c("chr_hg38",
               "MarkerName",
               "pos_hg38",
               "A1",
               "A2",
               "Freq",
               "LogOR",
               "StdErrLogOR",
               "P")

hg38_gwas <- hg38_gwas[, ..col_order]

names(hg38_gwas)[1] <- "CHR"

names(hg38_gwas)[3] <- "BP"

# read bim file with unique rsIDs
uniq_bim <-
    fread(
        "filter_data/Amygdala_unique_snps_bim/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus_dropBrains_maf01_hwe6_geno10_hg38_filtered_Amygdala_MDD_uniqueSNPs.bim"
    )

colnames(uniq_bim) <- c("CHR", "SNP", "dummy", "BP", "A1", "A2")

# check for overlap between gwas and unique bim file
uniq_bim$new_test <- paste0(uniq_bim$CHR, "_", uniq_bim$BP)

hg38_gwas$new_test <-
    paste0(gsub("chr", "" , hg38_gwas$CHR), "_", hg38_gwas$BP)

table(uniq_bim$new_test %in% hg38_gwas$new_test)

# > table(uniq_bim$new_test %in% hg38_gwas$new_test)
# 
# FALSE    TRUE 
# 4409489 6577690 
# > 4409489 / 657769

# merge in unique rsIDs
hg38_gwas <-
    merge(hg38_gwas[,-"SNP"], uniq_bim[, c("SNP", "new_test")], by = "new_test")

hg38_gwas$new_test <- NULL

col_order <- c(
    "CHR",
    "SNP",
    "BP",
    "A1",
    "A2",
    "Freq",
    "LogOR",
    "StdErrLogOR",
    "P")

hg38_gwas <- hg38_gwas[, ..col_order]

# hg38_gwas$effect <- log(hg38_gwas$OR)
hg38_gwas$Z <- hg38_gwas$LogOR / hg38_gwas$StdErrLogOR

pdf(file = "PGC_UKB_23andMe_depression_genome-wide_HIST.pdf", useDingbats = FALSE)

hist(hg38_gwas$effect, color = "gold")
hist(hg38_gwas$Z, color = "darkorange")

dev.off()

hg38_gwas_clean <-
    hg38_gwas[, c("SNP", "A1", "A2", "P", "Z")]
# colnames(hg38_gwas_clean)[4] <- "N"

write.table(
    hg38_gwas_clean,
    file = "PGC_UKB_23andMe_depression_genome-wide_CLEAN.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)


