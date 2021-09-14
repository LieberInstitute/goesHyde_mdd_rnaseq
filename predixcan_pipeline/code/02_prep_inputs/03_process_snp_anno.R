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
## Process genotype file
snp_gen <- snp

snp_gen$CHR <- snpMap$chr_hg38

snp_gen$CHR <- snp_gen$CHR %>% gsub("chr", "", .)

snp_gen$varID <- paste0(snpMap$chr_hg38, ":", snpMap$pos_hg38, ":", snpMap$COUNTED, ":", snpMap$ALT)
snp_gen$rsid <- snpMap$rsNumGuess

snp_gen <- snp_gen %>%
  relocate(rsid) %>%
  relocate(varID)

snp_gen <- snp_gen[!is.na(snp_gen$rsid),]

# table(snp_gen$varID %in% snp_anno$varID)
#
# FALSE   TRUE
# 154 410801

# head(snp_gen[!(snp_gen$varID  %in%  snp_anno$varID),1:5])
# varID        rsid 4463344439_R01C02 4463344373_R01C02
# 41391  NA:NA:T:C rs185317146                 2                 1
# 363431 NA:NA:T:C rs541293479                 0                 0
# 363452 NA:NA:A:T  rs61787357                 0                 0
# 369147 NA:NA:A:G  rs74813420                 0                 0
# 370016 NA:NA:T:G rs201128975                 1                 0
# 988980 NA:NA:T:C  rs62652648                 0                 0
# 4572348328_R01C01
# 41391                  1
# 363431                 0
# 363452                 0
# 369147                 0
# 370016                 0
# 988980                 0

split_snp_geno <- split(snp_gen, snp_gen$CHR)



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
