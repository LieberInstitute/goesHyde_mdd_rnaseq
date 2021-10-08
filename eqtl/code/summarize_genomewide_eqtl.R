library(SummarizedExperiment)
library(VariantAnnotation)
library(tidyverse)
library(jaffelab)
library(ggforce)
library(sessioninfo)
library(here)

#### genomewide results ####
## load tables
geneEqtl_tensor = map(c(amyg = "Amygdala", sacc = "sACC"),~read_csv(here("eqtl","data","tensorQTL_out",paste0("gene_",.x,"_cis_genomewide.csv"))))
geneEqtl_tensor$amyg

## matrixEQTL out
load(here("eqtl", "data", "matrixEQTL_out", "matrixEqtl_output_amyg_genomewide_gene_annotate.rda"), verbose = TRUE)
geneEqtl_matrix <- geneEqtl
table(geneEqtl_matrix$FDR < 0.01)
# FALSE     TRUE 
# 10602310  1574482 

geneEqtl_methods <- geneEqtl_tensor$amyg %>%
  select(snps = variant_id, gene = phenotype_id, statistic_tensor = slope, starts_with("pval")) %>%
  mutate(FDR_tensor = p.adjust(pval_perm, 'fdr')) %>%
  inner_join(as_tibble(geneEqtl_matrix) %>% rename(statistic_matrix = statistic, FDR_matrix = FDR))

geneEqtl_methods %>%
  ggplot(aes(FDR_matrix, FDR_tensor)) +
  geom_point()

with(geneEqtl_methods, table(FDR_matrix < 0.05, FDR_tensor < 0.05))

## list significant SNPs to make subset VCF
map(genomewide, ~.x %>% count(pval_perm < 0.01))
genomewide_signif <- do.call("rbind", map(genomewide, ~.x %>% filter(pval_perm < 0.05)))

signif_snps <- unique(genomewide_signif$variant_id)
length(signif_snps)
# 4,688 SNPs are significant
cat(signif_snps, file = here("eqtl","data","tensorQTL_out","significant_snps.txt"), sep = "\n")

## load subset VCF
signif_vcf <- readVcf(here("eqtl", "data", "signif_snps", "LIBD_Brain_merged_maf_005_topmed_051120_mdd_signif.vcf.gz"))
signif_vcf <- readGeno(here("eqtl", "data", "signif_snps", "LIBD_Brain_merged_maf_005_topmed_051120_mdd_signif.vcf.gz"), x = 'GT')
signif_info <- readInfo(here("eqtl", "data", "signif_snps", "LIBD_Brain_merged_maf_005_topmed_051120_mdd_signif.vcf.gz"), x = 'RS')
corner(signif_vcf)
head(signif_info)

