library(dplyr)
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

snp_anno <- snpMap[!is.na(snpMap$rsNumGuess),]

snp_anno <- snp_anno[, c("chr_hg38", "pos_hg38", "SNP", "COUNTED", "ALT")]

snp_anno$chr_hg38 <- snp_anno$chr_hg38 %>% gsub("chr", "", .)

snp_anno$snp_id_originalVCF <-
  paste0("snp_", gsub("chr", "", snp_anno$chr_hg38), "_", snp_anno$pos_hg38)

snp_anno$rsNumGuess <-
  snpMap[!is.na(snpMap$rsNumGuess),]$rsNumGuess

snp_anno$Num_alt_per_site <- nchar(snp_anno$ALT)

# colnames(snp_anno) <- c("Chr", "Pos", "VariantID", "Ref_b37", "Alt", "snp_id_originalVCF", "Num_alt_per_site")
colnames(snp_anno) <-
  c("chromosome",
    "pos",
    "varID",
    "ref_vcf",
    "alt_vcf",
    "snp_id_originalVCF",
    "rsid")

snp_anno <- snp_anno[snp_anno$chromosome != "X",]

# snp_anno$VariantID <- snp_anno$VariantID %>% gsub(pattern = "[:]", replacement = "_") %>% gsub(pattern = "chr", replacement = "")

snp_gen <- snp

snp_gen$varID <- snpMap$rsNumGuess

snp_gen <- snp_gen %>%
  relocate(varID)

snp_gen <- snp_gen[!is.na(snp_gen$varID),]

snp_gen$CHR <- snpMap[!is.na(snpMap$rsNumGuess),]$CHR

split_snp_geno <- split(snp_gen, snp_gen$CHR)

# > sapply(1:23, function(i) table(split_snp_geno[[i]]$CHR))
# 1    10    11    12    13    14    15    16    17    18    19     2    20
# 24707 17091 20370 15627 11830 10655  9276 12784  8675 17533 12203 29672  7674
# 21    22     3     4     5     6     7     8     9     X
# 4431  5019 23445 24640 21159 24902 21004 20506 67675    77

# for(i in 1:22){
#   split_snp_geno[[i]]$CHR <- NULL
# }

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