library("tidyverse")
library("here")
library("readxl")

## list all the SNPs of interest + signif to subset the VCF via plink

## SNPs of Intrest
twas_table <- read_excel(here("eqtl","data","plot_data","TWAS_sig_PP4_coloc08.xlsx")) %>%
  rename(variant_id_gwas = "chr_bp_ref_alt...6", variant_id_eqtl = "chr_bp_ref_alt...9")

snps_twas <- unique(c(twas_table$variant_id_gwas, twas_table$variant_id_eqtl))
length(snps_twas)
# [1] 49


#### Load the FDR01 output ####
fn_FDR01 <- list.files(here("eqtl", "data", "tensorQTL_FDR01"), recursive = TRUE, pattern = "*.csv", full.names = TRUE)
names(fn_FDR01) <- gsub("_FDR01.csv","", basename(fn_FDR01))

eqtl_FDR01 <- map(fn_FDR01, read.csv)
map_int(eqtl_FDR01, nrow)
# cis_gene_amyg         cis_gene_sacc independent_gene_amyg independent_gene_sacc     nominal_exon_amyg 
# 1326                  1367                  1805                  1880               5599118 
# nominal_exon_sacc     nominal_gene_amyg     nominal_gene_sacc      nominal_jxn_amyg      nominal_jxn_sacc 
# 5449203                625776                649396               7422175               7935775 
# nominal_tx_amyg       nominal_tx_sacc 
# 1365378               1444336 

sum(map_int(eqtl_FDR01, nrow))
# [1] 30497535


#### find all unique SNPs aka variant_id in the data ####
snps_FDR01 <- map(eqtl_FDR01, ~unique(.x$variant_id))
map_int(snps_FDR01, length)

snps_FDR01_unique <- unique(unlist(snps_FDR01))
length(snps_FDR01_unique)
# 2742184

all_snps <- unique(c(snps_FDR01_unique, snps_twas))
length(all_snps)
# [1] 2742188

## write to txt file
cat(all_snps, sep = "\n", file = here("eqtl", "data", "signif_snps",'significant_snps.txt'))
## run subset_signif_snps.sh to subset vcf
