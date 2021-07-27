load("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/goesHyde_bipolarMdd_Genotypes.rda")

snp_anno <- snpMap[, c("CHR", "POS", "SNP", "COUNTED", "ALT")]

snp_anno$snp_id_originalVCF <- paste0("snp_", snp_anno$CHR, "_", snp_anno$POS)

snp_anno$RSID_dbSNP137 <- snp_anno$rsNumGuess

snp_anno$Num_alt_per_site <- nchar(snp_anno$ALT)

colnames(snp_anno) <- c("Chr", "Pos", "VariantID", "Ref_b37", "Alt", "snp_id_originalVCF")

write.table(snp_anno, 
            file = here::here("predixcan_pipeline", "processed-data", "02_prep_inputs", "snp_annot_prep.txt"), 
            quote = FALSE,
            sep = "\t")