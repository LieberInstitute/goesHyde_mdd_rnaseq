library(tidyverse)
library(sessioninfo)
library(here)
library(miniparquet)

## Check out maf = 0 vs. maf = 0.05
eqtl <- parquet_read(here("eqtl", "data", "tensorQTL_out","genomewide_nominal","gene_Amygdala.cis_qtl_pairs.chr1.parquet")) %>%
  mutate(FDR = p.adjust(pval_nominal, "fdr"))
nrow(eqtl)
# [1] 8159237
summary(eqtl$af)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#     0.0     0.0     0.1     0.1     0.2     0.5  378453
table(eqtl$FDR < 0.05)
# FALSE    TRUE 
# 7241852   66965

## Results with maf filter
eqtl_maf <- parquet_read(here("eqtl", "data", "tensorQTL_out","genomewide_nominal_maf","gene_Amygdala.cis_qtl_pairs.chr1.parquet"))%>%
  mutate(FDR_maf = p.adjust(pval_nominal, "fdr"))
nrow(eqtl_maf)
# [1] 4554608
summary(eqtl_maf$af)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0500  0.1102  0.2111  0.2319  0.3417  0.515
nrow(eqtl_maf)/nrow(eqtl)
# [1] 0.5582149
table(eqtl_maf$FDR_maf < 0.05)
# FALSE    TRUE 
# 4504457   50151

eqtl_maf <- eqtl_maf %>% rename(pval_nominal_maf = pval_nominal, slope_maf = slope, slope_se_maf = slope_se)


eqtl_all <- left_join(eqtl_maf, eqtl, by = c("phenotype_id", "variant_id", "tss_distance")) %>%
  mutate(FDR_maf2 = p.adjust(pval_nominal, "fdr"))

head(eqtl_all)
nrow(eqtl_all)
# [1] 4554608
## compare values
eqtl_all %>% summary()

eqtl_all %>% filter(is.na(slope)) %>% head()

summary(eqtl_all$af.x - eqtl_all$af.y)
summary(eqtl_all$slope_maf - eqtl_all$slope)

## plots
slope_scatter <- ggplot(eqtl_all, aes(slope, slope_maf))+
  geom_point()

ggsave(slope_scatter, filename = here("eqtl", "plots", "test", "maf_slope_scatter.png"))

pval_scatter <- ggplot(eqtl_all, aes(pval_nominal, pval_nominal_maf))+
  geom_point()

ggsave(pval_scatter, filename = here("eqtl", "plots", "test", "maf_pval_scatter.png"))

FDR_scatter <- eqtl_all %>% 
  select(phenotype_id, variant_id, starts_with("FDR")) %>%
  pivot_longer(cols = c("FDR", "FDR_maf2"), names_to = 'FDR_no_maf', values_to = 'FDR') %>%
  mutate(`FDR no maf` = ifelse(grepl("maf2", FDR_no_maf),"Manual af filter FDR","All pair FDR")) %>%
  ggplot(aes(FDR_maf, FDR, color = `FDR no maf`))+
  geom_point()

ggsave(FDR_scatter, filename = here("eqtl", "plots", "test", "maf_FDR_scatter.png"))

eqtl_all %>% 
  count(FDR_maf < 0.05, 
        FDR_maf2 < 0.05)

#   FDR_maf < 0.05 FDR_maf2 < 0.05       n
# 1          FALSE           FALSE 4504453
# 2          FALSE            TRUE       4
# 3           TRUE            TRUE   50151
