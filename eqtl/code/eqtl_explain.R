library(SummarizedExperiment)
library(VariantAnnotation)
library(tidyverse)
library(jaffelab)
library(ggforce)
library(sessioninfo)
library(here)
library(miniparquet)
library(ggrepel)
# source(here("eqtl", "code", "utils.R"))
# load(here("data", "MDD_colors.Rdata"), verbose = TRUE)

#### load expression data ####
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)

rd <- as.data.frame(rowData(rse_gene)) %>% select(gencodeID, Symbol)

rr <- rowRanges(rse_gene) ## first 10 genes

g10 <- rr10$gencodeID


## read in AMYG eqtl data
nom <- read.csv(here("eqtl", "data", "tensorQTL_FDR01", "genomewide_nominal","nominal_gene_amyg_FDR01.csv"))
cis <- read.csv(here("eqtl", "data", "tensorQTL_out", "genomewide_cis","cis_gene_Amygdala.csv"))%>%
  rename(FDR = qval)
ind <- read.csv(here("eqtl", "data", "tensorQTL_out", "genomewide_independent","independent_gene_Amygdala.csv"), row.names = 1) %>%
  mutate(FDR = p.adjust(pval_perm, "fdr"))

nrow(nom)
# [1] 625776
nrow(cis)
# [1] 25085
nrow(ind)
# [1] 2225

## compare FDRs 
compare_FDR <- cis%>%
  select(phenotype_id, variant_id, cisFDR = FDR) %>%
  full_join(ind %>%
              select(phenotype_id, variant_id, indFDR = FDR)) %>%
  full_join(nom %>%
              select(phenotype_id, variant_id, nomFDR = FDR))

nrow(compare_FDR) #[1] [1] 638376

library(GGally)

FDR_scatter <- compare_FDR %>% 
  filter(!is.na(cisFDR), !is.na(indFDR), !is.na(nomFDR)) %>%
  ggpairs(columns= 3:5)

ggsave(FDR_scatter, filename = here("eqtl", "plots", "test", "FDR_scatter.png"))

FDR_scatter_all <- compare_FDR %>% 
  ggpairs(columns= 3:5)

ggsave(FDR_scatter_all, filename = here("eqtl", "plots", "test", "FDR_scatter_all.png"))


## find an example gene
nom_chr10 <- parquet_read(here("eqtl","data","tensorQTL_out", "genomewide_nominal","gene_Amygdala.cis_qtl_pairs.chr10.parquet"))%>%
  mutate(FDR = p.adjust(pval_nominal, "fdr")) ## this is just for demo, use full genome for read FDR correction!

nom_all_example <- nom_chr10 %>%
  filter(phenotype_id == "ENSG00000067057.16")

nrow(nom_all_example)
head(nom_all_example)


example <- compare_FDR %>%
  full_join(nom_all_example %>% select(phenotype_id, variant_id,tss_distance, nom1FDR = FDR)) %>%
  filter(phenotype_id == "ENSG00000067057.16") 

example_long <- example %>%
  pivot_longer(cols = ends_with("FDR"), 
               values_to = "FDR", 
               names_to = "eqtl")%>%
  separate(variant_id, into = c("chr", "bp","ref","alt"),remove = FALSE) %>%
  mutate(eqtl = gsub("FDR","",eqtl),
         bp = as.integer(bp), 
         eqtl2 = case_when(eqtl == "nom1"~"Nominal",
                           eqtl == "ind"~"Independent",
                           TRUE~eqtl)) %>%
  filter(!is.na(FDR), eqtl != "nom")

example_long %>% count(eqtl2, FDR < 0.01)

rr_example <- rr[rr$gencodeID == "ENSG00000067057.16",]

color_anno <- example_long %>% 
  filter(FDR < 0.01) %>%
  group_by(variant_id) %>%
  summarize(signif_eqtl = paste(gsub("1","",eqtl), collapse = " + "))

example_long$eqtl2 <- ordered(example_long$eqtl2, levels = c("Nominal", "cis", "Independent"))

example_long %>%
  count(eqtl2)
# eqtl2           n
# <ord>         <int>
# 1 Nominal      4637
# 2 cis             1
# 3 Independent     2

eqtl_example <- example_long %>%
  left_join(color_anno) %>%
  # replace_na(list(signif_eqtl = "none")) %>%
  ggplot(aes(x = bp, y = -log10(FDR), color = signif_eqtl)) +
  geom_point(alpha = 0.5) +
  geom_text_repel(aes(label= ifelse(signif_eqtl != "none", variant_id, "")), size = 3)+
  facet_wrap(~eqtl2)+
  geom_vline(xintercept = c(3066333, 3137712), linetype = 'dotted') +
  geom_hline(yintercept = -log10(0.01), color = "red", linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(eqtl_example, filename = here("eqtl", "plots", "test", "eqtl_example.png"), width = 10)
