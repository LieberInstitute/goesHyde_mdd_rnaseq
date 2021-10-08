library(SummarizedExperiment)
library(VariantAnnotation)
library(tidyverse)
library(jaffelab)
library(ggforce)
library(sessioninfo)
library(here)

source("utils.R")
load(here("data","MDD_colors.Rdata"), verbose = TRUE)

#### genomewide results ####
## load tables
geneEqtl_tensor = map(c(amyg = "Amygdala", sacc = "sACC"),~read.csv(here("eqtl","data","tensorQTL_out",paste0("gene_",.x,"_cis_genomewide.csv"))) %>%
                        mutate(FDR = p.adjust(pval_beta, 'fdr')) %>%
                        rename(gencodeID = phenotype_id))
head(geneEqtl_tensor$amyg)
summary(geneEqtl_tensor$amyg$tss_distance)

#### Compare Direct and Beta ####
## recommended on http://fastqtl.sourceforge.net/ "Checking that the Experiment Went Well"
p_val_check <- map(geneEqtl_tensor, ~ggplot(.x, aes(pval_perm, pval_beta)) + 
                     geom_point(size = 0.5, alpha = 0.5) +
                     geom_abline(color = "red"))

walk2(p_val_check, names(p_val_check), ~ggsave(.x, filename = here("eqtl", "plots", paste0("gene_genomewide_pval_check_", .y, ".png"))))

## top EQTLs
map(geneEqtl_tensor, ~.x %>% arrange(FDR) %>% head())
map_int(geneEqtl_tensor, ~.x %>% filter(FDR < 0.05) %>% nrow())
# amyg sacc 
# 2060 2092 
map_int(geneEqtl_tensor, ~.x %>% filter(FDR < 0.01) %>% nrow())

## matrixEQTL out
load(here("eqtl", "data", "matrixEQTL_out", "matrixEqtl_output_amyg_genomewide_gene_annotate.rda"), verbose = TRUE)
geneEqtl_matrix <- geneEqtl
table(geneEqtl_matrix$FDR < 0.01)
# FALSE     TRUE 
# 10602310  1574482 


geneEqtl_methods <- geneEqtl_tensor$amyg %>%
  select(snps = variant_id, gene = phenotype_id, statistic_tensor = slope, pval_beta) %>%
  mutate(FDR_tensor = p.adjust(pval_beta, 'fdr')) %>%
  inner_join(as_tibble(geneEqtl_matrix) %>% rename(statistic_matrix = statistic, FDR_matrix = FDR))

method_scater <- geneEqtl_methods %>%
  ggplot(aes(FDR_matrix, FDR_tensor)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_hline(yintercept = 0.05, color = "red") +
  geom_vline(xintercept = 0.05, color = "blue")

ggsave(method_scater, filename = here("eqtl", "plots","compare_methods_genomewide.png"))

with(geneEqtl_methods, table(FDR_matrix < 0.05, FDR_tensor < 0.05))

## list significant SNPs to make subset VCF ####
map(genomewide, ~.x %>% count(pval_perm < 0.01))
genomewide_signif <- do.call("rbind", map(genomewide, ~.x %>% filter(pval_perm < 0.05)))

signif_snps <- unique(genomewide_signif$variant_id)
length(signif_snps)
# 4,688 SNPs are significant
cat(signif_snps, file = here("eqtl","data","tensorQTL_out","significant_snps.txt"), sep = "\n")

#### load expression data ####
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rd <- as.data.frame(rowData(rse_gene)) %>% select(gencodeID, Symbol)

## build model
pd <- colData(rse_gene) %>% as.data.frame()
mod <- model.matrix(~PrimaryDx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, data = pd)
colnames(mod)

## get residual expression
gene_resid_expres_split <- map(splitit(rse_gene$BrainRegion), ~get_resid_expres(rse_gene[,.x], mod[.x,]))
map(gene_resid_expres_split, dim)

gene_resid_expres_split_long <- map(gene_resid_expres_split, ~.x %>%
                                      as.data.frame() %>%
                                      rownames_to_column("gencodeID") %>%
                                      pivot_longer(!gencodeID, names_to = "genoSample", values_to = "expression"))

## get phenotype data
pd <- as.data.frame(colData(rse_gene))
genoDx <- pd %>% select(genoSample, PrimaryDx) %>% unique()
dim(genoDx)

#### Load VCF ####
signif_vcf <- readVcf(here("eqtl", "data", "signif_snps", "LIBD_Brain_merged_maf_005_topmed_051120_mdd_signif.vcf.gz"))

rs_id <- info(signif_vcf ) %>%
  as.data.frame() %>%
  dplyr::select(RS) %>% 
  mutate(RS = paste0("rs", RS)) %>%
  rownames_to_column("variant_id")

head(rs_id)

geno_long <- as.data.frame(geno(signif_vcf)$GT) %>% 
  rownames_to_column("variant_id") %>%
  pivot_longer(!variant_id, names_to = "genoSample", values_to = "Genotype") %>%
  mutate(Genotype = case_when(Genotype == "0|0"~0,
                              Genotype == "1|1"~2,
                              TRUE~1),
         Genotype = as.factor(Genotype))

#### Filter and Annotate eQTL ####
geneEqtl_tensor_anno <- map(geneEqtl_tensor, function(e) {
 e_anno <- e %>% 
    filter(FDR < 0.01) %>%
    arrange(FDR) %>%
    left_join(rd) %>%
    left_join(rs_id) %>%
    mutate(eqtl = paste0(row_number(),". ", gencodeID, " - " ,variant_id),
           eqtl_anno = paste("Gene:",Symbol,"\nSNP:",RS, "\nFDR= ", scales::scientific(FDR, didgits = 3))) %>%
    slice(1:10)
 
 e_anno$eqtl <- factor(e_anno$eqtl)
 e_anno$eqtl <- fct_reorder(e_anno$eqtl, e_anno$FDR, min)

 return(e_anno)
 })


express_geno_cf <- map2(geneEqtl_tensor_anno, gene_resid_expres_split_long,
                        ~.x %>% 
                          left_join(geno_long, by = "variant_id") %>%
                          left_join(.y, by = c("gencodeID", "genoSample")) %>%
                          left_join(genoDx, by = "genoSample"))

map(express_geno_cf, head)

genomewide_box <- map2(express_geno_cf, geneEqtl_tensor_anno,
    ~ggplot(.x, aes(x = Genotype, y = expression)) +
    geom_boxplot(aes(fill = PrimaryDx)) + 
    facet_wrap(~eqtl, nrow = 2) +
    scale_fill_manual(values = mdd_Dx_colors) +
    geom_text(data = .y,
              aes(label = eqtl_anno),
              x = -Inf, y =Inf, vjust = "inward", hjust = "inward", size = 3, nudge_y = -0.1) +
    theme_bw()+ 
    coord_cartesian(clip = "off")
)

walk2(genomewide_box, names(genomewide_box), ~ggsave(.x, filename = here("eqtl", "plots",paste0("gene_genomewide_boxplot_",.y,".png")), width = 12))

