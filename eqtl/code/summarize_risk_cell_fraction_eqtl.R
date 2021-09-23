
library(tidyverse)
library(sessioninfo)
library(here)

tensorqtl_fn <- list.files(here("eqtl","data","tensorQTL_out","nominal_bpd_risk"), pattern = "*.csv", full.names = TRUE)

tensorqtl_out <- map_dfr(tensorqtl_fn, ~read_csv(.x) %>% mutate(fn = .x)) %>% 
  mutate(fn = gsub("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_out/nominal_bpd_risk/","",fn)) %>%
  separate(fn, into = c("feature", "region", "cell_type", NA))

summary(tensorqtl_out$tss_distance)

##846 Pairs each analysis
tensorqtl_out %>% group_by(region, cell_type) %>% count()

tensorqtl_adj <- tensorqtl_out %>%
  group_by(region, cell_type) %>%
  select(phenotype_id, variant_id, tss_distance, region, cell_type, pval_i, pval_gi) %>%
  mutate(FDR_gi_ct = p.adjust(pval_gi, "fdr"),
         p.bonf_gi_ct = p.adjust(pval_gi, "bonf")) %>%
  left_join(tensorqtl_out %>%
              group_by(region) %>%
              select(phenotype_id, variant_id, tss_distance, region, cell_type, pval_i, pval_gi) %>%
              mutate(FDR_gi_region = p.adjust(pval_gi, "fdr"),
                     p.bonf_gi_region = p.adjust(pval_gi, "bonf")))


tensorqtl_adj %>% filter(FDR_gi_ct < 0.05) %>% arrange(FDR_gi_ct)
tensorqtl_adj %>% filter(FDR_gi_ct < 0.05) %>% count() %>% arrange(-n)

sig_ct <- tensorqtl_adj  %>% summarise(FDR_ct_sig = sum(FDR_gi_ct < 0.05)) %>%
  left_join(tensorqtl_adj %>% summarise(bonf_ct_sig = sum(p.bonf_gi_ct < 0.05)))

sig_region <- tensorqtl_adj  %>% summarise(FDR_region_sig = sum(FDR_gi_region < 0.05)) %>%
  left_join(tensorqtl_adj %>% summarise(bonf_region_sig = sum(p.bonf_gi_region < 0.05)))


interaction_summary <- left_join(sig_ct, sig_region)

write_csv(interaction_summary, file = here("eqtl","data","summary","gene_BPD_risk_cell_fraction_interaction_summary.csv"))
