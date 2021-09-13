library(dplyr)
library(stringr)
# qrsh -l mem_free=50G,h_vmem=50G

setwd("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline")

load(
  here::here(
    "predixcan_pipeline",
    "processed-data",
    "02_prep_inputs",
    "goesHyde_bipolarMdd_Genotypes_PredictDB_NO-MDS.rda"
  )
)

# Remove all rows from snp annotation where there is no RSID
snp_anno <- snpMap[!is.na(snpMap$rsNumGuess),]

snp_anno <- snp_anno[, c("chr_hg38", "pos_hg38", "SNP", "COUNTED", "ALT", "rsNumGuess")]

# table(paste0(snp_anno$chr_hg38, ":", snp_anno$pos_hg38, ":", snp_anno$COUNTED, ":", snp_anno$ALT) == snp_anno$SNP)

snp_anno$SNP <- paste0(snp_anno$chr_hg38, ":", snp_anno$pos_hg38, ":", snp_anno$COUNTED, ":", snp_anno$ALT)

names(snp_anno)[3] <- "SNP_hg38"

snp_anno$chr_hg38 <- snp_anno$chr_hg38 %>% gsub("chr", "", .)

snp_anno$snp_id_originalVCF <-
  paste0("snp_", gsub("chr", "", snp_anno$chr_hg38), "_", snp_anno$pos_hg38)

snp_anno$Num_alt_per_site <- nchar(snp_anno$ALT)

# colnames(snp_anno) <- c("Chr", "Pos", "VariantID", "Ref_b37", "Alt", "snp_id_originalVCF", "Num_alt_per_site")
colnames(snp_anno) <-
  c("chromosome",
    "pos",
    "varID",
    "ref_vcf",
    "alt_vcf",
    "rsid",
    "snp_id_originalVCF",
    "Num_alt_per_site")

# TODO Should I be removing X?
snp_anno <- snp_anno[snp_anno$chromosome != "X",]

snp_anno <- snp_anno[!is.na(snp_anno),]

# snp_anno$VariantID <- snp_anno$VariantID %>% gsub(pattern = "[:]", replacement = "_") %>% gsub(pattern = "chr", replacement = "")

# Assuming that snp and snpMap are in the same order
snp_gen <- snp

snp_gen$CHR <- snpMap$chr_hg38

snp_anno$chr_hg38 <- snp_gen$CHR %>% gsub("chr", "", .)

snp_gen$varID <- paste0(snpMap$chr_hg38, ":", snpMap$pos_hg38, ":", snpMap$COUNTED, ":", snpMap$ALT)
snp_gen$rsid <- snpMap$rsNumGuess

snp_gen <- snp_gen %>%
  relocate(rsid) %>%
  relocate(varID)

snp_gen <- snp_gen[!is.na(snp_gen$rsid),]

# > table(snp_anno$varID %in% snp_gen$varID)
#
# FALSE    TRUE
# 2875607  410801
# > table(snp_gen$varID  %in%  snp_anno$varID)
#
# FALSE   TRUE
# 154 410801

split_snp_geno <- split(snp_gen, snp_gen$CHR)

# test1 <- sapply(split_snp_geno[[1]]$varID, str_split, ":")
#
# # TODO Why is this happening?
# table(sapply(test1, "[[", 1) %>% unname())

# > sapply(1:23, function(i) table(split_snp_geno[[i]]$CHR))
# 1    10    11    12    13    14    15    16    17    18    19     2    20
# 24707 17091 20370 15627 11830 10655  9276 12784  8675 17533 12203 29672  7674
# 21    22     3     4     5     6     7     8     9     X
# 4431  5019 23445 24640 21159 24902 21004 20506 67675    77

# str_split(split_snp_geno[[1]]$varID, ":")[[1]][1]


#
#  chr1 chr21  chr9    NA
# 24700     1     1     5

sapply(1:22,
       function (x)
         write.table(
           split_snp_geno[[x]],
           file = here::here(
             "predixcan_pipeline",
             "processed-data",
             "02_prep_inputs",
             "split_geno",
             paste0("split_snp_geno.chr",
                    split_snp_geno[[x]]$CHR[1], ".txt")
           ),
           sep = "\t",
           quote = FALSE,
           row.names = FALSE
         ))


write.table(
  snp_anno,
  file = here::here(
    "predixcan_pipeline",
    "processed-data",
    "02_prep_inputs",
    "snp_annot_prep.txt"
  ),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

# write.table(
#   snp_gen,
#   file = here::here(
#     "predixcan_pipeline",
#     "processed-data",
#     "02_prep_inputs",
#     "genotype.txt"
#   ),
#   quote = FALSE,
#   sep = "\t",
#   row.names = FALSE
# )
