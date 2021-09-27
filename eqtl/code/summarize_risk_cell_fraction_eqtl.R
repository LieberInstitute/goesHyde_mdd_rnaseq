library(SummarizedExperiment)
library(VariantAnnotation)
library(tidyverse)
library(jaffelab)
library(ggforce)
library(sessioninfo)
library(here)

## Load gene data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rd <- as.data.frame(rowData(rse_gene)) %>% select(gencodeID, Symbol)

## Load SNP data
risk_vcf <- readVcf(here("eqtl","data","risk_snps","LIBD_maf01_gwas_BPD.vcf.gz"))
rs_id <- info(risk_vcf) %>% as.data.frame() %>% dplyr::select(RS) %>% mutate(RS = paste0("rs", RS)) %>% rownames_to_column("variant_id")
head(rs_id)

geno_long <- as.data.frame(geno(risk_vcf)$GT) %>% 
  rownames_to_column("variant_id") %>%
  pivot_longer(!variant_id, names_to = "genoSample", values_to = "Genotype") %>%
  mutate(Genotype = case_when(Genotype == "0|0"~"0",
                              Genotype == "1|1"~"2",
                              TRUE~"1"))

## Load tensorqtl outputs
# amyg_astro <- read.csv(here("eqtl","data","tensorQTL_out","nominal_bpd_risk","gene_Amygdala_Astro.csv"))
# amyg_astro_mm <- read.csv(here("eqtl","data","tensorQTL_out","nominal_bpd_risk","gene_Amygdala_Astro_mm.csv"))
# summary(amyg_astro$pval_gi - amyg_astro_mm$pval_gi)

tensorqtl_fn <- list.files(here("eqtl","data","tensorQTL_out","nominal_bpd_risk"), pattern = "*mm.csv", full.names = TRUE)

tensorqtl_out <- map_dfr(tensorqtl_fn, ~read.csv(.x) %>% mutate(fn = .x)) %>% 
  mutate(fn = gsub("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_out/nominal_bpd_risk/","",fn)) %>%
  separate(fn, into = c("feature", "region", "cell_type", NA)) %>%
  rename(gencodeID = phenotype_id)

summary(tensorqtl_out$tss_distance)

##846 Pairs each analysis
tensorqtl_out %>% group_by(region, cell_type) %>% count()

tensorqtl_adj <- tensorqtl_out %>%
  group_by(region, cell_type)  %>%
  left_join(rd, by = "gencodeID") %>%
  left_join(rs_id, by = "variant_id") %>%
  dplyr::select(BrainRegion = region, gencodeID, Symbol, variant_id, RS, cell_type, tss_distance, pval_i, pval_gi) %>%
  mutate(FDR_gi = p.adjust(pval_gi, "fdr"),
         p.bonf_gi = p.adjust(pval_gi, "bonf"))

tensorqtl_adj %>% filter(FDR_gi < 0.05) %>% arrange(FDR_gi)
tensorqtl_adj %>% filter(FDR_gi < 0.05) %>% ungroup() %>% dplyr::count(BrainRegion)
tensorqtl_adj %>% filter(FDR_gi < 0.05) %>% count() %>% arrange(-n)
tensorqtl_adj %>% filter(FDR_gi < 0.05) %>% ungroup() %>% count(BrainRegion, gencodeID, variant_id) %>% arrange(-n)

write_csv(tensorqtl_adj, file = here("eqtl","data","summary","gene_BPD_risk_cell_fraction_interaction.csv"))

## Summarise number of significant
interaction_summary <- tensorqtl_adj  %>% summarise(FDR_ct_sig = sum(FDR_gi < 0.05)) %>%
  left_join(tensorqtl_adj %>% summarise(bonf_ct_sig = sum(p.bonf_gi < 0.05)))

write_csv(interaction_summary, file = here("eqtl","data","summary","gene_BPD_risk_cell_fraction_interaction_summary.csv"))

#### Prep Data for Plotting ####

get_resid_expres <- function(rse, mod, tpm = FALSE){
  mod <- mod[colnames(rse),]
  if(tpm) rpkm <- assays(rse)$tpm else rpkm <- recount::getRPKM(rse, "Length")
  ## residualize expression
  exprs <- log2(rpkm + 1)
  exprs <- cleaningY(exprs, mod, P = 1)
  colnames(exprs) <- rse$genoSample
  return(exprs)
}

## build model
pd <- colData(rse_gene) %>% as.data.frame()
mod <- model.matrix(~PrimaryDx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, data = pd)
colnames(mod)

gene_resid_expres <- get_resid_expres(rse_gene, mod)
corner(gene_resid_expres)

gene_resid_long <- gene_resid_expres %>% 
  as.data.frame() %>%
  rownames_to_column("gencodeID") %>%
  pivot_longer(!gencodeID, names_to = "genoSample", values_to = "expression")

## Cell fractions
cell_fractions <- pd %>% 
  dplyr::select(BrainRegion, genoSample, Astro, Endo, Macro, Micro, Mural, Oligo, OPC, Tcell, Excit, Inhib) %>%
  pivot_longer(!c(BrainRegion, genoSample), names_to = "cell_type", values_to = "cell_fraction")

## Merge Data
tensorqtl_adj_anno <- tensorqtl_adj %>%
  filter(FDR_gi < 0.05) %>%
  ungroup() %>% 
  dplyr::group_by(BrainRegion) %>% 
  arrange(FDR_gi) %>%
  mutate(eqtl = paste0(row_number(),". ", gencodeID, " - " ,variant_id),
         eqtl_anno = paste("Gene:",Symbol,"\nSNP:",RS, "\nFDR= ", scales::scientific(FDR_gi, didgits = 3)))

tensorqtl_adj_anno$eqtl <- factor(tensorqtl_adj_anno$eqtl)
tensorqtl_adj_anno$eqtl <- fct_reorder(tensorqtl_adj_anno$eqtl, tensorqtl_adj_anno$FDR_gi, min)

levels(tensorqtl_adj_anno$eqtl)
tensorqtl_adj_anno %>% count(eqtl) 

express_geno_cf <- tensorqtl_adj_anno %>% 
  left_join(geno_long, by = "variant_id") %>%
  left_join(cell_fractions, by = c("BrainRegion", "cell_type", "genoSample")) %>%
  left_join(gene_resid_long, by = c("gencodeID", "genoSample"))

express_geno_cf %>% count(cell_type, eqtl)

#### PLOT ####
tensorqtl_adj_anno %>% count(BrainRegion)

walk(c("sACC","Amygdala"), function(r){
  express_cf_sactter <- express_geno_cf %>%
    filter(BrainRegion == r) %>%
    ggplot(aes(x = cell_fraction, y = expression)) +
    geom_point(aes(color = Genotype), alpha = 0.5, size = .5) +
    geom_smooth(aes(color = Genotype, fill = Genotype), method = lm)+
    geom_text(aes(label = eqtl_anno), x = 0, y =Inf, vjust = "inward", hjust = "inward", size = 3)
  
  required_n_pages <- n_pages(express_cf_sactter +  facet_wrap_paginate(~eqtl + cell_type, scales = "free", nrow = 2, ncol =2))
  message(required_n_pages)
  
  
  pdf(here("eqtl","plots",paste0("eqtl_risk_BPD_cell_fraction_",r,".pdf")), width = 10)
  for(i in 1:required_n_pages){
    print(express_cf_sactter +  facet_wrap_paginate(~eqtl + cell_type, scales = "free", nrow = 2, ncol =2, page = i))
  }
  dev.off()
})

