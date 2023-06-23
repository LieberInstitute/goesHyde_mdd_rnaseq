library("SummarizedExperiment")
library("VariantAnnotation")
library("tidyverse")
library("jaffelab")
# library(ggforce)
library("sessioninfo")
library("here")
# library("xlsx")

source(here("eqtl", "code", "utils.R"))
load(here("data", "MDD_colors.Rdata"), verbose = TRUE)

# current_signif_snps <- readLines(here("eqtl", "data", "signif_snps",'significant_snps.txt'))
# length(current_signif_snps)

# twas_table <- read.xlsx(here("eqtl","data","plot_data","TWAS_sig_PP4_coloc08.xlsx"), sheetIndex = "Sheet1", header=TRUE)

## read.xlsx thorws warning?
# april_plots <- read.xlsx(here("eqtl","data","plot_data","Genes_for_boxplot_042622.xlsx"), sheetIndex = "eqtl", header=TRUE)

april_plots <- readxl::read_excel(here("eqtl","data","plot_data","Genes_for_boxplot_042622.xlsx"), sheet = "eqtl")
# all(april_plots$variant_id %in% current_signif_snps)
# [1] FALSE


#### load expression data ####
## move code to prep plot eqtl
load(here("eqtl", "data", "plot_data", "residual_expres.Rdata"), verbose = TRUE)
resid_expres_split <- resid_expres_split$gene
names(resid_expres_split)

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
