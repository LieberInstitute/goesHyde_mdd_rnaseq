library("SummarizedExperiment")
library("VariantAnnotation")
library("tidyverse")
library("jaffelab")
# library("ggforce")
library("sessioninfo")
library("here")

source(here("eqtl", "code", "utils.R"))
load(here("data", "MDD_colors.Rdata"), verbose = TRUE)

#### load expression data ####
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rd <- as.data.frame(rowData(rse_gene)) %>% select(gencodeID, Symbol)

# features <- c("gene", "exon", "jxn", "tx")
features <- c("gene")
names(features) <- features
regions <- c(amyg = "Amygdala", sacc = "sACC")

#### cis genomewide results ####
## gene only for test
nominal_files <- list.files(here("eqtl", "data", "tensorQTL_FDR05", "genomewide_nominal"),
    pattern = "gene", full.names = TRUE
)
names(nominal_files) <- gsub("_FDR05.csv", "", basename(nominal_files))

cisEqtl_tensor <- map(nominal_files, read.csv)

map(cisEqtl_tensor, head)

#### matrixEQTL out ####
load(here("eqtl", "data", "matrixEQTL_out", "matrixEqtl_output_amyg_genomewide_gene_annotate.rda"), verbose = TRUE)
eqtl_matrix <- geneEqtl %>%
    as.data.frame() %>%
    select(variant_id = snps, phenotype_id = gene, slope_ME = statistic, pval_ME = pvalue, FDR_ME = FDR)

table(eqtl_matrix$FDR_ME < 0.01)
# FALSE     TRUE
# 10602310  1574482
table(eqtl_matrix$FDR_ME < 0.05)

summary(eqtl_matrix$pval)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.005343 0.030638 0.036567 0.063498 0.100000 


#### Compare Gene Amyg
eqtl_methods <- cisEqtl_tensor$gene_amyg %>%
    inner_join(eqtl_matrix, by = c("phenotype_id", "variant_id"))


eqtl_methods %>% count(FDR < 0.01, FDR_ME< 0.01)
# FDR < 0.01 FDR_ME < 0.01      n
# 1      FALSE         FALSE  33667
# 2      FALSE          TRUE 105567
# 3       TRUE         FALSE  10560
# 4       TRUE          TRUE 319015

eqtl_methods <- eqtl_methods %>% 
    mutate(signif = case_when(FDR < 0.01 & FDR_ME < 0.01 ~ "Both",
                              FDR < 0.01 ~ "tensor",
                              FDR_ME < 0.01 ~ "matrix",
                              TRUE~"FDR > 0.01"
                              ))

eqtl_methods %>% count(signif)

nrow(eqtl_methods)
# [1] 468809

slope_scater <- eqtl_methods %>%
    ggplot(aes(x = slope_ME, y = slope, color = signif)) +
    geom_point(size = 0.5, alpha = 0.5) +
    labs(x = "slope MatrixEQTL", y = "slope tensorQTL") +
    facet_wrap(~signif)

ggsave(slope_scater, filename = here("eqtl", "plots", "compare_methods_slope.png"))


# ## list significant SNPs to make subset VCF ####
# map(genomewide, ~ .x %>% count(pval_perm < 0.01))
# genomewide_signif <- do.call("rbind", map(genomewide, ~ .x %>% filter(pval_perm < 0.05)))
# 
# signif_snps <- unique(genomewide_signif$variant_id)
# length(signif_snps)
# # 4,688 SNPs are significant
# cat(signif_snps, file = here("eqtl", "data", "tensorQTL_out", "significant_snps.txt"), sep = "\n")
# 
# 
# ## build model
# pd <- colData(rse_gene) %>% as.data.frame()
# mod <- model.matrix(~ PrimaryDx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, data = pd)
# colnames(mod)
# 
# ## get residual expression
# gene_resid_expres_split <- map(splitit(rse_gene$BrainRegion), ~ get_resid_expres(rse_gene[, .x], mod[.x, ]))
# map(gene_resid_expres_split, dim)
# 
# gene_resid_expres_split_long <- map(gene_resid_expres_split, ~ .x %>%
#     as.data.frame() %>%
#     rownames_to_column("gencodeID") %>%
#     pivot_longer(!gencodeID, names_to = "genoSample", values_to = "expression"))
# 
# ## get phenotype data
# pd <- as.data.frame(colData(rse_gene))
# genoDx <- pd %>%
#     select(genoSample, PrimaryDx) %>%
#     unique()
# dim(genoDx)
# 
# #### Load VCF ####
# signif_vcf <- readVcf(here("eqtl", "data", "signif_snps", "LIBD_Brain_merged_maf_005_topmed_051120_mdd_signif.vcf.gz"))
# 
# rs_id <- info(signif_vcf) %>%
#     as.data.frame() %>%
#     dplyr::select(RS) %>%
#     mutate(RS = paste0("rs", RS)) %>%
#     rownames_to_column("variant_id")
# 
# head(rs_id)
# 
# geno_long <- as.data.frame(geno(signif_vcf)$GT) %>%
#     rownames_to_column("variant_id") %>%
#     pivot_longer(!variant_id, names_to = "genoSample", values_to = "Genotype") %>%
#     mutate(
#         Genotype = case_when(
#             Genotype == "0|0" ~ 0,
#             Genotype == "1|1" ~ 2,
#             TRUE ~ 1
#         ),
#         Genotype = as.factor(Genotype)
#     )
# 
# #### Filter and Annotate eQTL ####
# geneEqtl_tensor_anno <- map(geneEqtl_tensor, function(e) {
#     e_anno <- e %>%
#         filter(FDR < 0.01) %>%
#         arrange(FDR) %>%
#         left_join(rd) %>%
#         left_join(rs_id) %>%
#         mutate(
#             eqtl = paste0(row_number(), ". ", gencodeID, " - ", variant_id),
#             eqtl_anno = paste("Gene:", Symbol, "\nSNP:", RS, "\nFDR= ", scales::scientific(FDR, didgits = 3))
#         ) %>%
#         slice(1:10)
# 
#     e_anno$eqtl <- factor(e_anno$eqtl)
#     e_anno$eqtl <- fct_reorder(e_anno$eqtl, e_anno$FDR, min)
# 
#     return(e_anno)
# })
# 
# 
# express_geno_cf <- map2(
#     geneEqtl_tensor_anno, gene_resid_expres_split_long,
#     ~ .x %>%
#         left_join(geno_long, by = "variant_id") %>%
#         left_join(.y, by = c("gencodeID", "genoSample")) %>%
#         left_join(genoDx, by = "genoSample")
# )
# 
# map(express_geno_cf, head)
# 
# genomewide_box <- map2(
#     express_geno_cf, geneEqtl_tensor_anno,
#     ~ ggplot(.x, aes(x = Genotype, y = expression)) +
#         geom_boxplot(aes(fill = PrimaryDx)) +
#         facet_wrap(~eqtl, nrow = 2) +
#         scale_fill_manual(values = mdd_Dx_colors) +
#         geom_text(
#             data = .y,
#             aes(label = eqtl_anno),
#             x = -Inf, y = Inf, vjust = "inward", hjust = "inward", size = 3, nudge_y = -0.1
#         ) +
#         theme_bw() +
#         coord_cartesian(clip = "off")
# )
# 
# walk2(genomewide_box, names(genomewide_box), ~ ggsave(.x, filename = here("eqtl", "plots", paste0("gene_genomewide_boxplot_", .y, ".png")), width = 12))

# sgejobs::job_single("compare_genomewide_eqtl", memory = "25G",create_shell = TRUE, command = "Rscript compare_genomewide_eqtl.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
