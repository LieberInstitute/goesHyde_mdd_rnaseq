library(SummarizedExperiment)
library(VariantAnnotation)
library(tidyverse)
library(jaffelab)
library(ggforce)
library(sessioninfo)
library(here)

source(utils.R)

#### Get Gene Residual Expression ####
## Load gene data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rd <- as.data.frame(rowData(rse_gene)) %>% select(gencodeID, Symbol)

## build model
pd <- colData(rse_gene) %>% as.data.frame()
mod <- model.matrix(~PrimaryDx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, data = pd)
colnames(mod)

## get residual expression by region
gene_resid_expres_split <- map(splitit(rse_gene$BrainRegion), ~get_resid_expres(rse_gene[,.x], mod[.x,]))
gene_resid_expres <- do.call("cbind", gene_resid_expres_split)

corner(gene_resid_expres)

gene_resid_long <- gene_resid_expres %>% 
  as.data.frame() %>%
  rownames_to_column("gencodeID") %>%
  pivot_longer(!gencodeID, names_to = "genoSample", values_to = "expression")

#### BPD Risk Data ####
## Load SNP data
risk_vcf <- readVcf(here("eqtl","data","risk_snps","LIBD_maf01_gwas_BPD.vcf.gz"))
rs_id <- info(risk_vcf) %>%
  as.data.frame() %>%
  dplyr::select(RS) %>% 
  mutate(RS = paste0("rs", RS)) %>%
  rownames_to_column("variant_id")

## Risk SNP details
risk_snps <- read_csv(here("eqtl","data","risk_snps","pgc3-supptable2.csv"))
summary(risk_snps$`OR (A1 allele)`)

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
## good!

geno_long <- as.data.frame(geno(risk_vcf)$GT) %>% 
  rownames_to_column("variant_id") %>%
  pivot_longer(!variant_id, names_to = "genoSample", values_to = "Genotype") %>%
  mutate(Genotype = case_when(Genotype == "0|0"~0,
                              Genotype == "1|1"~2,
                              TRUE~1),
         Genotype = as.factor(Genotype))

## Load tensorqtl outputs
tensorqtl_fn <- list.files(here("eqtl","data","tensorQTL_out","nominal_bpd_risk"), pattern = "*.csv", full.names = TRUE)

tensorqtl_out <- map_dfr(tensorqtl_fn, ~read.csv(.x) %>% mutate(fn = .x)) %>% 
  mutate(fn = gsub("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_out/nominal_bpd_risk/","",fn)) %>%
  separate(fn, into = c("feature", "region", "cell_type", NA)) %>%
  rename(gencodeID = phenotype_id)

summary(tensorqtl_out$tss_distance)
## 1Mbp window

##846 Pairs each analysis
tensorqtl_out %>% group_by(region, cell_type) %>% count()
length(unique(tensorqtl_out$variant_id))
# [1] 59

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

## Cell fractions
cell_fractions <- pd %>% 
  dplyr::select(BrainRegion, genoSample, Astro, Endo, Macro, Micro, Mural, Oligo, OPC, Tcell, Excit, Inhib) %>%
  pivot_longer(!c(BrainRegion, genoSample), names_to = "cell_type", values_to = "cell_fraction")

## Merge Data
tensorqtl_adj_anno <- tensorqtl_adj %>%
  filter(FDR_gi < 0.05) %>%
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

express_geno_anno <- tensorqtl_adj_anno %>% 
  select(BrainRegion, variant_id, eqtl, eqtl_anno, cell_type) %>%
  left_join(risk_snps_detials) %>% 
  mutate(eqtl_anno = paste0(eqtl_anno, "\n", risk_anno))

#### PLOT ####
tensorqtl_adj_anno %>% count(BrainRegion)

walk(c("Amygdala", "sACC"), function(r){
  message(r)
  express_cf_sactter <- express_geno_cf %>%
    filter(BrainRegion == r) %>%
    ggplot(aes(x = cell_fraction, y = expression)) +
    geom_point(aes(color = Genotype), alpha = 0.5, size = .5) +
    geom_smooth(aes(color = Genotype, fill = Genotype), method = lm) 
  
  required_n_pages <- n_pages(express_cf_sactter +  facet_wrap_paginate(~eqtl + cell_type, scales = "free", nrow = 2, ncol =2))
  message(required_n_pages)
  
  pdf(here("eqtl","plots",paste0("eqtl_risk_BPD_cell_fraction_",r,".pdf")), width = 10)
  for(i in 1:required_n_pages){
    print(express_cf_sactter +  
            facet_wrap_paginate(~eqtl + cell_type, scales = "free", nrow = 2, ncol =2, page = i) 
          +
            geom_text(data = filter(express_geno_anno, BrainRegion == r),
                      aes(label = eqtl_anno), x = -Inf, y =Inf, vjust = "inward", hjust = "inward", size = 3)
          )
  }
  dev.off()
})

