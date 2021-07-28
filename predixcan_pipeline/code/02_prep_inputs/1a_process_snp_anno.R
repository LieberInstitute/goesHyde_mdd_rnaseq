library(dplyr)
# qrsh -l mem_free=50G,h_vmem=50G

load("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/goesHyde_bipolarMdd_Genotypes.rda")

snp_anno <- snpMap[, c("CHR", "POS", "SNP", "COUNTED", "ALT")]

# snp_anno$snp_id_originalVCF <- paste0("snp_", snp_anno$CHR, "_", snp_anno$POS)

snp_anno$RSID_dbSNP137 <- snp_anno$rsNumGuess

# snp_anno$Num_alt_per_site <- nchar(snp_anno$ALT)

# colnames(snp_anno) <- c("Chr", "Pos", "VariantID", "Ref_b37", "Alt", "snp_id_originalVCF", "Num_alt_per_site")
colnames(snp_anno) <- c("chromosome", "Pos", "VariantID", "ref_vcf", "alt_vcf")

# snp_anno$VariantID <- snp_anno$VariantID %>% gsub(pattern = "[:]", replacement = "_") %>% gsub(pattern = "chr", replacement = "")

snp_gen <- snp

snp_gen$VariantID <- snp_anno$VariantID

snp_gen <- snp_gen %>% 
  relocate(VariantID)

write.table(snp_anno, 
            file = here::here("predixcan_pipeline", "processed-data", "02_prep_inputs", "snp_annot_prep.txt"), 
            quote = FALSE,
            sep = "\t")

write.table(snp_gen, 
            file = here::here("predixcan_pipeline", "processed-data", "02_prep_inputs", "genotype.txt"), 
            quote = FALSE,
            sep = "\t")