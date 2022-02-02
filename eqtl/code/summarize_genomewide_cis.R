library(SummarizedExperiment)
library(VariantAnnotation)
library(tidyverse)
library(jaffelab)
library(ggforce)
library(sessioninfo)
library(here)

source(here("eqtl", "code", "utils.R"))
load(here("data", "MDD_colors.Rdata"), verbose = TRUE)

summarize_eqtl <- function(eqtl_df) {
    eqtl_stats <- list(
        n_snp_feature_pairs = nrow(eqtl_df),
        n_snps = length(unique(eqtl_df$variant_id)),
        n_features = length(unique(eqtl_df$phenotype_id))
    )
    return(eqtl_stats)
}

#### load expression data ####
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
# load(here("exprs_cutoff", "rse_exon.Rdata"), verbose = TRUE)
# load(here("exprs_cutoff", "rse_jxn.Rdata"), verbose = TRUE)
# load(here("exprs_cutoff", "rse_tx.Rdata"), verbose = TRUE)

rd <- as.data.frame(rowData(rse_gene)) %>% select(gencodeID, Symbol)

features <- c("gene", "exon", "jxn", "tx")
names(features) <- features
regions <- c(amyg = "Amygdala", sacc = "sACC")

#### genomewide results ####
## load tables
message("Loading tensorQTL cis results")

eqtl_tensor <- map(
    c(amyg = "Amygdala", sacc = "sACC"),
    ~ read.csv(here("eqtl", "data", "tensorQTL_out", "genomewide_cis", paste0("gene_", .x, ".csv"))) %>%
        rename(
            gencodeID = phenotype_id,
            FDR = qval
        ) %>%
        separate(variant_id, into = c("chr", NA, NA, NA), sep = ":", remove = FALSE)
)


head(eqtl_tensor$amyg)
summary(eqtl_tensor$amyg$tss_distance)
summary(eqtl_tensor$amyg$FDR)

map_int(eqtl_tensor, nrow)
# amyg  sacc
# 25085 25085

map(eqtl_tensor, ~ .x %>% count(FDR < 0.05))
order <- map2(eqtl_tensor, c("sACC", "Amygdala"), ~ .x %>%
    filter(FDR < 0.01) %>%
    count(chr) %>%
    mutate(region = .y))

# order2 <- do.call("rbind", order) %>%
#     arrange(n)
#
# cat(paste0(order2$region, "_", order2$chr), file = here("eqtl","code","region_chr.txt"), sep = "\n")
# # $amyg
# FDR < 0.01     n
# 1      FALSE   321
# 2       TRUE 24764
#
# $sacc
# FDR < 0.01     n
# 1      FALSE   263
# 2       TRUE 24822

## Remember one result per gene
map_int(eqtl_tensor, ~ length(unique(.x$gencodeID)))
# amyg  sacc
# 25085 25085

map_int(eqtl_tensor, ~ length(unique(.x$variant_id)))
# amyg  sacc
# 16339 16027

#### How many are risk SNPs? ####
mdd_snps <- read.delim(here("eqtl", "data", "risk_snps", "PGC_MDD_genome-wide_significant_Jan2022.txt")) %>%
    filter(!is.na(bp)) %>%
    rename(variant_id = chr_bp_ref_alt)

mdd_snps %>% count(chr)

cis_risk <- map(eqtl_tensor, ~ .x %>% inner_join(mdd_snps))
map_int(cis_risk, nrow)

cis_risk2 <- do.call("rbind", cis_risk) %>%
    rownames_to_column("region") %>%
    mutate(region = gsub("\\.[0-9]+", "", region))

cis_risk2 %>%
    count(variant_id) %>%
    filter(n == 2)
cis_risk2 %>%
    count(gencodeID) %>%
    filter(n == 2)



#### Prep data for plotting ####
load(here("eqtl", "data", "plot_data", "residual_expres.Rdata"), verbose = TRUE)

resid_expres_split <- resid_expres_split$gene

resid_expres_split_long <- map(resid_expres_split, ~ .x %>%
    as.data.frame() %>%
    rownames_to_column("gencodeID") %>%
    pivot_longer(!gencodeID, names_to = "genoSample", values_to = "expression"))
map_int(resid_expres_split_long, nrow)

## get phenotype data
pd <- as.data.frame(colData(rse_gene))
genoDx <- pd %>%
    select(genoSample, PrimaryDx) %>%
    unique()
head(genoDx)

#### Load VCF ####
signif_vcf <- readVcf(here("eqtl", "data", "signif_snps", "LIBD_Brain_merged_maf_005_topmed_051120_mdd_signif.vcf.gz"))
dim(signif_vcf)

## connect variant id to rs number
rs_id <- info(signif_vcf) %>%
    as.data.frame() %>%
    dplyr::select(RS) %>%
    mutate(RS = paste0("rs", RS)) %>%
    rownames_to_column("variant_id")

map_int(eqtl_tensor, ~ length(unique(.x$variant_id)))
map_int(eqtl_tensor, ~ sum(unique(.x$variant_id) %in% rs_id$variant_id))

## get genotype for samples in long format
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
eqtl_tensor_top$gene$amyg %>%
    count(gencodeID) %>%
    arrange(n) %>%
    left_join(rd)

## annotation for boxplots - order eqtl pairs by FDR
eqtl_tensor_anno <- map_depth(eqtl_tensor_top, 2, function(e) {
    e_anno <- e %>%
        arrange(FDR) %>%
        left_join(rd) %>% ## add gene symbol
        left_join(rs_id) %>% ## add rs #
        mutate(
            eqtl = paste0(str_pad(row_number(), 4), ". ", gencodeID, " - ", variant_id),
            eqtl_anno = paste("\n Gene:", Symbol, "\n SNP:", RS, "\n FDR= ", scales::scientific(FDR, didgits = 3))
        ) %>%
        slice(c(1:50, 201:250, 401:450, 601:650, 801:850, 1001:1050))

    e_anno$eqtl <- factor(e_anno$eqtl)
    e_anno$eqtl <- fct_reorder(e_anno$eqtl, e_anno$FDR, min)

    return(e_anno)
})
map_depth(eqtl_tensor_anno, 2, head)
map_depth(eqtl_tensor_anno, 2, nrow)

## merge genotype and expression data for boxplots
express_geno <- map2(
    eqtl_tensor_anno, resid_expres_split_long, function(eqtl, expres) {
        map2(eqtl, expres, ~ .x %>%
            inner_join(geno_long, by = "variant_id") %>%
            inner_join(.y, by = c("gencodeID", "genoSample")) %>%
            left_join(genoDx, by = "genoSample"))
    }
)

pwalk(
    list(expres = express_geno, anno = eqtl_tensor_anno, feat_name = names(express_geno)), function(expres, anno, feat_name) {
        pwalk(
            list(expres2 = expres, anno2 = anno, region_name = names(expres)), function(expres2, anno2, region_name) {
                label <- paste0(feat_name, "_", region_name)
                message(label)

                eqtl_box <- ggplot(expres2, aes(x = Genotype, y = expression)) +
                    geom_boxplot(aes(fill = PrimaryDx)) +
                    scale_fill_manual(values = mdd_Dx_colors) +
                    geom_text(
                        data = anno2,
                        aes(label = eqtl_anno),
                        x = -Inf, y = Inf, vjust = "inward", hjust = "inward", size = 3
                    ) +
                    theme_bw() +
                    coord_cartesian(clip = "off")

                required_n_pages <- n_pages(eqtl_box + facet_wrap_paginate(~eqtl, scales = "free", nrow = 2, ncol = 2))
                message(paste("Plotting", required_n_pages, "pages"))

                pdf(here("eqtl", "plots", paste0("eqtl_genomewide_", label, ".pdf")), width = 10)
                for (i in 1:required_n_pages) {
                    print(eqtl_box +
                        facet_wrap_paginate(~eqtl, scales = "free", nrow = 2, ncol = 2, page = i))
                }
                dev.off()
            }
        )
    }
)


# sgejobs::job_single('summarize_genomewide_eqtl', create_shell = TRUE, memory = '100G', command = "Rscript summarize_genomewide_eqtl.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
