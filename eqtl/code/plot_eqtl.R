library(SummarizedExperiment)
library(VariantAnnotation)
library(tidyverse)
library(jaffelab)
# library(ggforce)
library(sessioninfo)
library(here)

source(here("eqtl", "code", "utils.R"))
load(here("data", "MDD_colors.Rdata"), verbose = TRUE)

#### load expression data ####
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)

rd <- as.data.frame(rowData(rse_gene)) %>% select(gencodeID, Symbol)
pd <- as.data.frame(colData(rse_gene))

## build model
mod <- model.matrix(~ PrimaryDx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, data = pd)
colnames(mod)

## get residual expression by region
gene_resid_expres_split <- map(splitit(rse_gene$BrainRegion), ~ get_resid_expres(rse_gene[, .x], mod[.x, ]))

gene_resid_long <- map2(gene_resid_expres_split, names(gene_resid_expres_split),
                        ~.x %>%
                          as.data.frame() %>%
                          rownames_to_column("gencodeID") %>%
                          pivot_longer(!gencodeID, names_to = "genoSample", values_to = "expression") %>%
                          mutate(BrainRegion = .y))

gene_resid_long <- do.call("rbind", gene_resid_long)

#### load VCF ####
# sftp://jhpce01.jhsph.edu/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/signif_snps/LIBD_Brain_merged_maf_005_topmed_051120_mdd_signif.vcf.gz
signif_vcf <- readVcf(here("eqtl", "data", "signif_snps","LIBD_Brain_merged_maf_005_topmed_051120_mdd_signif.vcf.gz"))
dim(signif_vcf)
info(signif_vcf)

rs_id <- info(signif_vcf) %>%
  as.data.frame() %>%
  dplyr::select(RS) %>%
  mutate(RS = paste0("rs", RS)) %>%
  rownames_to_column("variant_id")

geno_long <- as.data.frame(geno(signif_vcf)$GT) %>%
  rownames_to_column("variant_id") %>%
  pivot_longer(!variant_id, names_to = "genoSample", values_to = "Genotype") %>%
  mutate(
    Genotype = case_when(
      Genotype == "0|0" ~ 0,
      Genotype == "1|1" ~ 2,
      TRUE ~ 1
    ),
    Genotype = as.factor(Genotype)
  )


#### genomewide results ####
## load tables
regions <- c(amyg = "Amygdala", sacc = "sACC")

eqtl_out <- map2(regions, names(regions),
                ~ read.csv(here("eqtl", "data", "tensorQTL_FDR01", "genomewide_independent", paste0("independent_gene_", .y, "_FDR01.csv")),
                           row.names = 1) %>%
                  mutate(BrainRegion = .x)
)

map_int(eqtl_out, nrow)

eqtl_ind <- do.call("rbind", eqtl_out)%>%
  as_tibble() %>%
  group_by(BrainRegion) %>%
  arrange(FDR) %>%
  rename(gencodeID = phenotype_id) %>%
  left_join(rd) %>%
  left_join(rs_id) %>%
  mutate(
    eqtl = as.factor(paste0(row_number(), ". ", gencodeID, " - ", variant_id)),
    eqtl2 = paste0("Gene: ", gencodeID, "\nSNP: ", variant_id),
    eqtl_anno = paste("Gene:", Symbol, "\nSNP:", RS, "\nFDR= ", scales::scientific(FDR, didgits = 3))
  )

eqtl_ind$eqtl <- fct_reorder(eqtl_ind$eqtl, eqtl_ind$FDR, min)
head(levels(eqtl_ind$eqtl))


eqtl_ind %>%
  select(BrainRegion, eqtl, FDR, eqtl_anno) %>%
  group_by(BrainRegion) %>%
  slice(1:10)

## check snps on vcf 
ind_signif_snps <- unique(eqtl_ind$variant_id)
length(ind_signif_snps)

## need to remake the vcf
table(ind_signif_snps %in% rownames(signif_vcf))
# FALSE  TRUE 
# 2508   103 

## for now just use 10 w/ SNPs in VCF
eqtl_ind2 <- eqtl_ind %>%
  filter(!is.na(RS)) %>%
  slice(1:10)

eqtl_anno <- eqtl_ind2 %>%
  select(eqtl, BrainRegion, eqtl_anno)

express_geno <- eqtl_ind2 %>%
  select(gencodeID, variant_id, FDR, BrainRegion, Symbol, RS, eqtl, eqtl2, eqtl_anno) %>%
  left_join(gene_resid_long) %>%
  left_join(geno_long)

express_geno

eqtl_box <- express_geno %>%
  filter(BrainRegion == "Amygdala") %>%
  ggplot(aes(x = Genotype, y = expression, fill = Genotype)) +
  geom_boxplot()+
  facet_wrap(~eqtl2, scales = "free")+
  # scale_fill_manual() +
  theme_bw()

ggsave(eqtl_box, filename = here("eqtl","plots","test","ind_box_test.png"), width = 10)
