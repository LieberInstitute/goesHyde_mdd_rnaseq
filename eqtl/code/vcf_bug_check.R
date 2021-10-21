library(SummarizedExperiment)
library(VariantAnnotation)
library(tidyverse)
library(here)

#### BPD Risk Data ####
## Load SNP data
risk_vcf <- readVcf(here("eqtl","data","risk_snps","LIBD_maf01_gwas_BPD_bug.vcf.gz"))
rs_id <- info(risk_vcf) %>%
  as.data.frame() %>%
  dplyr::select(RS) %>% 
  mutate(RS = paste0("rs", RS)) %>%
  rownames_to_column("variant_id")

## Risk SNP details
risk_snps <- read_csv(here("eqtl","data","risk_snps","pgc3-supptable2.csv"))

risk_snps_detials <- risk_snps %>% 
  separate(`A1/A2`, into = c("A1", "A2"), sep = "/") %>%
  select(RS = SNP, A1, A2, OR = `OR (A1 allele)`) %>%
  left_join(rs_id) %>%
  mutate(risk_snp = ifelse(OR > 1, "ref", "alt"),
         risk_anno = paste0("OR = ", OR, " Risk Allel:", risk_snp))

#missing 1 SNP -n we expect this!
table(risk_snps$SNP %in% rs_id$RS)
## Extra SNPs?
table(rs_id$RS %in% risk_snps$SNP)
rs_id[!rs_id$RS %in% risk_snps$SNP,]
## How did these SNPs end up in the VCF?
#              variant_id          RS
# 13   chr6:25759784:CT:C  rs67712859 * 
# 25 chr12:80069577:TAA:T rs112481528 *
# 27     chr2:4832797:C:T rs115694470
# 59  chr13:65469547:A:AT  rs35306829

good_rs <- intersect(rs_id$RS, risk_snps$SNP)
length(good_rs)

good_ids <- rs_id$variant_id[rs_id$RS %in% risk_snps$SNP]
risk_vcf_good <- risk_vcf[good_ids,]
dim(risk_vcf_good)

writeVcf(risk_vcf_good ,here("eqtl","data","risk_snps","LIBD_maf01_gwas_BPD.vcf"))
# module load htslib
# bgzip

#### MDD risk vcf ####

## table has GRCh38 annotations get SNPids to subset
mdd_snps <- read.delim(here("eqtl", "data", "risk_snps", "PGC_depression_genome-wide_significant_makers.txt"))
head(mdd_snps)
# chr       bp markername a1 a2   freq   logor stderrlogor         p
# 1   1 17795514  rs4141983  T  C 0.6740  0.0264      0.0046 9.692e-09
# 2   1 49062731 rs17105472  A  C 0.0960 -0.0409      0.0073 2.355e-08
# 3   1 49088804   rs589436  T  C 0.0960 -0.0409      0.0073 2.333e-08
# 4   1 49097909   rs655065  T  C 0.0963 -0.0412      0.0073 1.708e-08
# 5   1 49115334   rs590503  T  C 0.9037  0.0416      0.0073 1.321e-08
# 6   1 49209604   rs354155  C  G 0.0923 -0.0449      0.0075 1.751e-09

mdd_snps <- mdd_snps %>%
  mutate(snpID = paste0("chr",chr,":" ,bp, ":",a1, ":",a2))

cat(mdd_snps$snpID, file = here("eqtl", "data", "risk_snps","MDD_risk_snps.txt"), sep = "\n")


