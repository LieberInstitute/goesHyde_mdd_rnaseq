library(SummarizedExperiment)
library(VariantAnnotation)
library(tidyverse)
library(jaffelab)
library(ggforce)
library(sessioninfo)
library(here)

source(here("eqtl", "code", "utils.R"))
load(here("data", "MDD_colors.Rdata"), verbose = TRUE)

#### load expression data ####
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rd <- as.data.frame(rowData(rse_gene)) %>% select(gencodeID, Symbol)

#### genomewide results ####
## load tables
geneEqtl_tensor <- map(c(amyg = "Amygdala", sacc = "sACC"), ~ read.csv(here("eqtl", "data", "tensorQTL_out", paste0("gene_", .x, "_cis_genomewide.csv"))) %>%
    mutate(FDR = p.adjust(pval_beta, "fdr")) %>%
    rename(gencodeID = phenotype_id))
head(geneEqtl_tensor$amyg)
summary(geneEqtl_tensor$amyg$tss_distance)

#### Compare Direct and Beta ####
## recommended on http://fastqtl.sourceforge.net/ on in "Permutation pass in Cis" page: "Checking that the Experiment Went Well"
p_val_check <- map(geneEqtl_tensor, ~ ggplot(.x, aes(pval_perm, pval_beta)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_abline(color = "red"))

walk2(p_val_check, names(p_val_check), ~ ggsave(.x, filename = here("eqtl", "plots", paste0("gene_genomewide_pval_check_", .y, ".png"))))

## top EQTLs
map(geneEqtl_tensor, ~ .x %>%
    arrange(FDR) %>%
    head())
map_int(geneEqtl_tensor, ~ .x %>%
    filter(FDR < 0.05) %>%
    nrow())
# amyg sacc
# 2060 2092
map_int(geneEqtl_tensor, ~ .x %>%
    filter(FDR < 0.01) %>%
    nrow())
# amyg sacc
# 1322 1350

#### map_nominal results####
nominal_dir <- here("eqtl", "data", "tensorQTL_out", "gene_Amygdala_cis_genomewide_nominal")
parquet_files <- list.files(nominal_dir, full.names = TRUE)

dim(geneEqtl_nominal)
# [1] 91717646        10
length(unique(geneEqtl_nominal$variant_id))
# [1] 9736432
table(geneEqtl_nominal$FDR_nominal < 0.05)
# FALSE     TRUE
# 81048381   946823

## save out results < 0.01
geneEqtl_nominal_filter <- geneEqtl_nominal %>% filter(pval_nominal < 0.01)
dim(geneEqtl_nominal_filter)
# [1] 2428478      10
write_csv(geneEqtl_nominal_filter, file = paste0(nominal_dir, ".csv"))

#### matrixEQTL out ####
load(here("eqtl", "data", "matrixEQTL_out", "matrixEqtl_output_amyg_genomewide_gene_annotate.rda"), verbose = TRUE)
geneEqtl_matrix <- geneEqtl
table(geneEqtl_matrix$FDR < 0.01)
# FALSE     TRUE
# 10602310  1574482

## N snps
length(unique(geneEqtl_matrix$snps))
# [1] 4826322
length(unique(geneEqtl_tensor$amyg$variant_id))
# [1] 15767
length(intersect(geneEqtl_matrix$snps, geneEqtl_tensor$amyg$variant_id))
# [1] 9079

## N genes
nrow(rse_gene)
# [1] 25212
length(unique(geneEqtl_matrix$gene))
# [1] 24470
length(unique(geneEqtl_tensor$amyg$gencodeID))
# [1] 24241
length(intersect(geneEqtl_matrix$gene, geneEqtl_tensor$amyg$gencodeID))
# [1] 24039

## p-vals
summary(geneEqtl_tensor$amyg$pval_beta)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0000  0.1374  0.4250  0.4293  0.6982  0.9982
summary(geneEqtl_matrix$pvalue)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000000 0.005343 0.030638 0.036567 0.063498 0.100000
summary(geneEqtl_matrix$FDR)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0000  0.1322  0.3789  0.3323  0.5235  0.6184

## match up
geneEqtl_methods <- geneEqtl_tensor$amyg %>%
    select(snps = variant_id, gene = gencodeID, statistic_tensor = slope, pval_beta) %>%
    mutate(FDR_tensor = p.adjust(pval_beta, "fdr")) %>%
    inner_join(as_tibble(geneEqtl_matrix) %>%
        rename(statistic_matrix = statistic, FDR_matrix = FDR))

nrow(geneEqtl_methods)
# [1] 5371

method_scater <- geneEqtl_methods %>%
    ggplot(aes(FDR_matrix, FDR_tensor)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_hline(yintercept = 0.05, color = "red") +
    geom_vline(xintercept = 0.05, color = "blue")

ggsave(method_scater, filename = here("eqtl", "plots", "compare_methods_genomewide.png"))

with(geneEqtl_methods, table(FDR_matrix < 0.05, FDR_tensor < 0.05))
#       FALSE TRUE
# FALSE  2382    0
# TRUE   2196  793

## matchup w/ nominal
geneEqtl_methods_nominal <- geneEqtl_nominal %>%
    select(snps = variant_id, gene = phenotype_id, slope_tensor = slope, pval_nominal, FDR_nominal) %>%
    inner_join(as_tibble(geneEqtl_matrix) %>%
        rename(statistic_matrix = statistic, FDR_matrix = FDR, pval_matrix = pvalue)) %>%
    mutate(
        log10_FDR_matrix = -log10(FDR_matrix),
        log10_FDR_nominal = -log10(FDR_nominal)
    )

dim(geneEqtl_methods_nominal)
# [1] 8539512      13

method_density_nom <- geneEqtl_methods_nominal %>%
    ggplot(aes(log10_FDR_matrix, log10_FDR_nominal)) +
    geom_point(size = 0.01, alpha = 0.1) +
    geom_abline(color = "red") +
    theme_bw() +
    labs(x = "-log10(FDR matrixEQTL)", y = "-log10(FDR tensorQTL nominal)")

ggsave(method_density_nom, filename = here("eqtl", "plots", "compare_methods_genomewide_nom.png"))

with(geneEqtl_methods_nominal, table("nom" = FDR_nominal < 0.01, "matrix" = FDR_matrix < 0.01))
#         FALSE    TRUE
# FALSE 7015864 1176873
# TRUE    15681  331094

## list significant SNPs to make subset VCF ####
map(genomewide, ~ .x %>% count(pval_perm < 0.01))
genomewide_signif <- do.call("rbind", map(genomewide, ~ .x %>% filter(pval_perm < 0.05)))

signif_snps <- unique(genomewide_signif$variant_id)
length(signif_snps)
# 4,688 SNPs are significant
cat(signif_snps, file = here("eqtl", "data", "tensorQTL_out", "significant_snps.txt"), sep = "\n")


## build model
pd <- colData(rse_gene) %>% as.data.frame()
mod <- model.matrix(~ PrimaryDx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, data = pd)
colnames(mod)

## get residual expression
gene_resid_expres_split <- map(splitit(rse_gene$BrainRegion), ~ get_resid_expres(rse_gene[, .x], mod[.x, ]))
map(gene_resid_expres_split, dim)

gene_resid_expres_split_long <- map(gene_resid_expres_split, ~ .x %>%
    as.data.frame() %>%
    rownames_to_column("gencodeID") %>%
    pivot_longer(!gencodeID, names_to = "genoSample", values_to = "expression"))

## get phenotype data
pd <- as.data.frame(colData(rse_gene))
genoDx <- pd %>%
    select(genoSample, PrimaryDx) %>%
    unique()
dim(genoDx)

#### Load VCF ####
signif_vcf <- readVcf(here("eqtl", "data", "signif_snps", "LIBD_Brain_merged_maf_005_topmed_051120_mdd_signif.vcf.gz"))

rs_id <- info(signif_vcf) %>%
    as.data.frame() %>%
    dplyr::select(RS) %>%
    mutate(RS = paste0("rs", RS)) %>%
    rownames_to_column("variant_id")

head(rs_id)

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

#### Filter and Annotate eQTL ####
geneEqtl_tensor_anno <- map(geneEqtl_tensor, function(e) {
    e_anno <- e %>%
        filter(FDR < 0.01) %>%
        arrange(FDR) %>%
        left_join(rd) %>%
        left_join(rs_id) %>%
        mutate(
            eqtl = paste0(row_number(), ". ", gencodeID, " - ", variant_id),
            eqtl_anno = paste("Gene:", Symbol, "\nSNP:", RS, "\nFDR= ", scales::scientific(FDR, didgits = 3))
        ) %>%
        slice(1:10)

    e_anno$eqtl <- factor(e_anno$eqtl)
    e_anno$eqtl <- fct_reorder(e_anno$eqtl, e_anno$FDR, min)

    return(e_anno)
})


express_geno_cf <- map2(
    geneEqtl_tensor_anno, gene_resid_expres_split_long,
    ~ .x %>%
        left_join(geno_long, by = "variant_id") %>%
        left_join(.y, by = c("gencodeID", "genoSample")) %>%
        left_join(genoDx, by = "genoSample")
)

map(express_geno_cf, head)

genomewide_box <- map2(
    express_geno_cf, geneEqtl_tensor_anno,
    ~ ggplot(.x, aes(x = Genotype, y = expression)) +
        geom_boxplot(aes(fill = PrimaryDx)) +
        facet_wrap(~eqtl, nrow = 2) +
        scale_fill_manual(values = mdd_Dx_colors) +
        geom_text(
            data = .y,
            aes(label = eqtl_anno),
            x = -Inf, y = Inf, vjust = "inward", hjust = "inward", size = 3, nudge_y = -0.1
        ) +
        theme_bw() +
        coord_cartesian(clip = "off")
)

walk2(genomewide_box, names(genomewide_box), ~ ggsave(.x, filename = here("eqtl", "plots", paste0("gene_genomewide_boxplot_", .y, ".png")), width = 12))
