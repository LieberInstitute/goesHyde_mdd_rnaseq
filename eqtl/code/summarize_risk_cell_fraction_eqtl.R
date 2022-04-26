library("SummarizedExperiment")
library("VariantAnnotation")
library("tidyverse")
library("jaffelab")
library("ggforce")
library("sessioninfo")
library("here")

# source("utils.R")

#### Read in csv files ####
## loop over dx
dx <- c("mdd", "bpd")
names(dx) <- dx

regions <- c(amyg = "Amygdala", sacc = "sACC")

## get data from csv files
combos <- cross2(dx, names(regions))
names(combos) <- map_chr(cross2(dx, names(regions)),  ~paste0(.x[[1]],"_",.x[[2]]))

eqtl_out_risk <- map(combos, function(f){
    fn <- here("eqtl", "data", "tensorQTL_FDR05", paste0("risk_nominal_", f[[1]]),
               paste0("risk_gene_",f[[2]],"_",f[[1]],"_FDRall.csv"))
    
    eqtl_out <- read.csv(fn) %>%
        mutate(region = f[[2]],
               Dx = f[[1]],
               .before = "cell_type")
    
    return(eqtl_out)
})

map(eqtl_out_risk, head)

map(eqtl_out_risk, ~.x %>% 
        group_by(cell_type) %>% 
        summarize(n_pairs = n(),
                  n_FDR05 = sum(FDR < 0.05)))


map(eqtl_out_risk, ~.x %>% 
        group_by(cell_type, phenotype_id) %>% 
        filter(FDR < 0.05) %>%
        summarize(n_SNP = n(),
                  ))

## Combine in to one df
eqtl_risk_all <- do.call("rbind", eqtl_out_risk)
rownames(eqtl_risk_all) <- NULL

risk_summary <- eqtl_risk_all %>%
    filter(FDR < 0.05) %>%
    group_by(region, Dx, cell_type) %>%
    summarize(n_pairs = n(),
              n_genes = length(unique(phenotype_id)),
              n_SNPs = length(unique(variant_id))) %>%
    mutate(anno = paste("pairs:", n_pairs, "\ngenes:",n_genes,"\nSNPs:",n_SNPs))

risk_summary %>%
    select(region, Dx, cell_type, n_pairs) %>%
    pivot_wider(names_from = "cell_type", values_from = "n_pairs")

risk_summary %>%
    select(region, Dx, cell_type, n_genes) %>%
    pivot_wider(names_from = "cell_type", values_from = "n_genes")

risk_summary %>%
    select(region, Dx, cell_type, n_SNPs) %>%
    pivot_wider(names_from = "cell_type", values_from = "n_SNPs")

tile_plot <- risk_summary %>%
    ggplot(aes(x = region, y = cell_type, fill = n_pairs))+
    geom_tile() +
    geom_text(aes(label = anno), color = "gray") +
    facet_wrap(~Dx) +
    theme_bw() +
    labs(title = "Risk SNP cell type interaction")

ggsave(tile_plot, file = here("eqtl", "plots","eqtl_risk_n-pairs_tile.png"))

#### Get Gene Residual Expression ####
## Load gene data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rd <- as.data.frame(rowData(rse_gene)) %>% select(gencodeID, Symbol)
pd <- as.data.frame(colData(rse_gene))
## Cell fractions
cell_fractions <- pd %>%
    dplyr::select(BrainRegion, genoSample, Astro, Endo, Macro, Micro, Mural, Oligo, OPC, Tcell, Excit, Inhib) %>%
    pivot_longer(!c(BrainRegion, genoSample), names_to = "cell_type", values_to = "cell_fraction")

## build model
pd <- colData(rse_gene) %>% as.data.frame()
mod <- model.matrix(~ PrimaryDx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, data = pd)
colnames(mod)

## get residual expression by region
gene_resid_expres_split <- map(splitit(rse_gene$BrainRegion), ~ get_resid_expres(rse_gene[, .x], mod[.x, ]))
gene_resid_expres <- do.call("cbind", gene_resid_expres_split)

corner(gene_resid_expres)

gene_resid_long <- gene_resid_expres %>%
    as.data.frame() %>%
    rownames_to_column("gencodeID") %>%
    pivot_longer(!gencodeID, names_to = "genoSample", values_to = "expression")

#### BPD Risk Data ####
## Load SNP data
risk_vcf_bpd <- readVcf(here("eqtl", "data", "risk_snps", "LIBD_maf01_gwas_BPD.vcf.gz"))
rs_id_bpd <- info(risk_vcf_bpd) %>%
    as.data.frame() %>%
    dplyr::select(RS) %>%
    mutate(RS = paste0("rs", RS)) %>%
    rownames_to_column("variant_id")

## Risk SNP details
risk_snps_bpd <- read_csv(here("eqtl", "data", "risk_snps", "pgc3-supptable2.csv"))
summary(risk_snps_bpd$`OR (A1 allele)`)

risk_snps_bpd_detials <- risk_snps_bpd %>%
    separate(`A1/A2`, into = c("A1", "A2"), sep = "/") %>%
    select(RS = SNP, A1, A2, OR = `OR (A1 allele)`) %>%
    left_join(rs_id_bpd) %>%
    mutate(
        risk_snp = ifelse(OR > 1, "ref", "alt"),
        risk_anno = paste0("OR = ", OR, " Risk Allel:", risk_snp)
    )

# missing 1 SNP -n we expect this!
table(risk_snps_bpd$SNP %in% rs_id_bpd$RS)
## Extra SNPs?
table(rs_id_bpd$RS %in% risk_snps_bpd$SNP)
## good!

geno_long_bpd <- as.data.frame(geno(risk_vcf_bpd)$GT) %>%
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

## Load tensorqtl outputs
tensorqtl_fn_bpd <- list.files(here("eqtl", "data", "tensorqtl_out_bpd", "nominal_bpd_risk"), pattern = "*.csv", full.names = TRUE)

tensorqtl_out_bpd <- map_dfr(tensorqtl_fn_bpd, ~ read.csv(.x) %>% mutate(fn = .x)) %>%
    mutate(fn = gsub("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/tensorqtl_out_bpd/nominal_bpd_risk/", "", fn)) %>%
    separate(fn, into = c("feature", "region", "cell_type", NA)) %>%
    rename(gencodeID = phenotype_id)

summary(tensorqtl_out_bpd$tss_distance)
## 1Mbp window

## 846 Pairs each analysis
tensorqtl_out_bpd %>%
    group_by(region, cell_type) %>%
    count()
length(unique(tensorqtl_out_bpd$variant_id))
# [1] 59

tensorqtl_adj_bpd <- tensorqtl_out_bpd %>%
    group_by(region, cell_type) %>%
    left_join(rd, by = "gencodeID") %>%
    left_join(rs_id_bpd, by = "variant_id") %>%
    dplyr::select(BrainRegion = region, gencodeID, Symbol, variant_id, RS, cell_type, tss_distance, pval_i, pval_gi) %>%
    mutate(
        FDR_gi = p.adjust(pval_gi, "fdr"),
        p.bonf_gi = p.adjust(pval_gi, "bonf")
    )

tensorqtl_adj_bpd %>%
    filter(FDR_gi < 0.05) %>%
    arrange(FDR_gi)
tensorqtl_adj_bpd %>%
    filter(FDR_gi < 0.05) %>%
    ungroup() %>%
    dplyr::count(BrainRegion)
tensorqtl_adj_bpd %>%
    filter(FDR_gi < 0.05) %>%
    count() %>%
    arrange(-n)
tensorqtl_adj_bpd %>%
    filter(FDR_gi < 0.05) %>%
    ungroup() %>%
    count(BrainRegion, gencodeID, variant_id) %>%
    arrange(-n)

write_csv(tensorqtl_adj_bpd, file = here("eqtl", "data", "summary", "gene_BPD_risk_cell_fraction_interaction.csv"))

## Summarise number of significant
interaction_summary_bpd <- tensorqtl_adj_bpd %>%
    summarise(FDR_ct_sig = sum(FDR_gi < 0.05)) %>%
    left_join(tensorqtl_adj_bpd %>% summarise(bonf_ct_sig = sum(p.bonf_gi < 0.05)))

write_csv(interaction_summary_bpd, file = here("eqtl", "data", "summary", "gene_BPD_risk_cell_fraction_interaction_summary.csv"))

#### Prep BPD Data for Plotting ####
## Merge Data
tensorqtl_adj_bpd_anno <- tensorqtl_adj_bpd %>%
    filter(FDR_gi < 0.05) %>%
    dplyr::group_by(BrainRegion) %>%
    arrange(FDR_gi) %>%
    mutate(
        eqtl = paste0(row_number(), ". ", gencodeID, " - ", variant_id),
        eqtl_anno = paste("Gene:", Symbol, "\nSNP:", RS, "\nFDR= ", scales::scientific(FDR_gi, didgits = 3))
    )

tensorqtl_adj_bpd_anno$eqtl <- factor(tensorqtl_adj_bpd_anno$eqtl)
tensorqtl_adj_bpd_anno$eqtl <- fct_reorder(tensorqtl_adj_bpd_anno$eqtl, tensorqtl_adj_bpd_anno$FDR_gi, min)

levels(tensorqtl_adj_bpd_anno$eqtl)
tensorqtl_adj_bpd_anno %>% count(eqtl)

express_geno_bpd_cf <- tensorqtl_adj_bpd_anno %>%
    left_join(geno_long_bpd, by = "variant_id") %>%
    left_join(cell_fractions, by = c("BrainRegion", "cell_type", "genoSample")) %>%
    left_join(gene_resid_long, by = c("gencodeID", "genoSample"))

express_geno_anno <- tensorqtl_adj_bpd_anno %>%
    select(BrainRegion, variant_id, eqtl, eqtl_anno, cell_type) %>%
    left_join(risk_snps_bpd_detials) %>%
    mutate(eqtl_anno = paste0(eqtl_anno, "\n", risk_anno))

#### PLOT ####
tensorqtl_adj_bpd_anno %>% count(BrainRegion)

walk(c("Amygdala", "sACC"), function(r) {
    message(r)
    express_cf_sactter <- express_geno_cf %>%
        filter(BrainRegion == r) %>%
        ggplot(aes(x = cell_fraction, y = expression)) +
        geom_point(aes(color = Genotype), alpha = 0.5, size = .5) +
        geom_smooth(aes(color = Genotype, fill = Genotype), method = lm)

    required_n_pages <- n_pages(express_cf_sactter + facet_wrap_paginate(~ eqtl + cell_type, scales = "free", nrow = 2, ncol = 2))
    message(required_n_pages)

    pdf(here("eqtl", "plots", paste0("eqtl_risk_BPD_cell_fraction_", r, ".pdf")), width = 10)
    for (i in 1:required_n_pages) {
        print(express_cf_sactter +
            facet_wrap_paginate(~ eqtl + cell_type, scales = "free", nrow = 2, ncol = 2, page = i)
            +
            geom_text(
                data = filter(express_geno_anno, BrainRegion == r),
                aes(label = eqtl_anno), x = -Inf, y = Inf, vjust = "inward", hjust = "inward", size = 3
            ))
    }
    dev.off()
})

#### MDD risk snps ####
tensorqtl_fn_mdd <- list.files(here("eqtl", "data", "tensorQTL_out", "nominal_mdd_risk"), pattern = "*.csv", full.names = TRUE)

tensorqtl_out_mdd <- map_dfr(tensorqtl_fn_mdd, ~ read.csv(.x) %>% mutate(fn = .x)) %>%
    mutate(fn = basename(fn)) %>%
    separate(fn, into = c("feature", "region", "cell_type", NA)) %>%
    rename(gencodeID = phenotype_id)

dim(tensorqtl_out_mdd)
# [1] 24600    18
head(tensorqtl_out_mdd)
summary(tensorqtl_out_mdd$tss_distance)

# ## Test eign results
# tensorqtl_out_mdd_eigen <- read.csv(tensorqtl_fn_mdd[1]) %>%
#     select(gencodeID = phenotype_id, variant_id, pval_gi_eigen = pval_gi, b_gi_eigen = b_gi)
#
# dim(tensorqtl_out_mdd_eigen) #155
#
# eign_test <- tensorqtl_out_mdd %>%
#     filter(region == "Amygdala", cell_type == "Astro") %>%
#     inner_join(tensorqtl_out_mdd_eigen)
#
# dim(eign_test)
# ## seem to get 1:1 results, cool
# eign_test %>% ggplot(aes(b_gi, b_gi_eigen)) + geom_point()


## Load SNP data
risk_vcf_mdd <- readVcf(here("eqtl", "data", "risk_snps", "LIBD_maf01_gwas_MDD.vcf.gz"))
rs_id_mdd <- info(risk_vcf_mdd) %>%
    as.data.frame() %>%
    dplyr::select(RS) %>%
    mutate(RS = paste0("rs", RS)) %>%
    rownames_to_column("variant_id")

risk_snps_mdd <- read_delim(here("eqtl", "data", "risk_snps", "PGC_MDD_genome-wide_significant_Jan2022.txt")) %>%
    filter(!is.na(bp))

table(risk_snps_mdd$rsid %in% rs_id_mdd$RS)
# TRUE 
# 9650

risk_snps_mdd_detials <- risk_snps_mdd %>%
    select(RS = rsid, LogOR) %>%
    mutate(OR = exp(LogOR)) %>%
    inner_join(rs_id_mdd) %>%
    mutate(
        risk_snp = ifelse(OR > 1, "ref", "alt"),
        risk_anno = paste0("OR = ", round(OR, 3), " Risk Allel:", risk_snp)
    )

geno_long_mdd <- as.data.frame(geno(risk_vcf_mdd)$GT) %>%
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

## summarize
tensorqtl_adj_mdd <- tensorqtl_out_mdd %>%
    group_by(region, cell_type) %>%
    left_join(rd, by = "gencodeID") %>%
    left_join(rs_id_mdd, by = "variant_id") %>%
    dplyr::select(BrainRegion = region, gencodeID, Symbol, variant_id, RS, cell_type, tss_distance, pval_i, pval_gi) %>%
    mutate(
        FDR_gi = p.adjust(pval_gi, "fdr"),
        p.bonf_gi = p.adjust(pval_gi, "bonf")
    )

tensorqtl_adj_mdd %>%
    filter(FDR_gi < 0.05) %>%
    arrange(FDR_gi)

tensorqtl_adj_mdd %>%
    filter(FDR_gi < 0.05) %>%
    ungroup() %>%
    dplyr::count(BrainRegion)

# # A tibble: 2 × 2
# BrainRegion     n
# <chr>       <int>
# 1 Amygdala       52
# 2 sACC          118

tensorqtl_adj_mdd %>%
    filter(FDR_gi < 0.05) %>%
    count() %>%
    arrange(-n)

tensorqtl_adj_mdd %>%
    filter(FDR_gi < 0.05) %>%
    ungroup() %>%
    count(BrainRegion, gencodeID, variant_id) %>%
    arrange(-n)

write_csv(tensorqtl_adj_mdd, file = here("eqtl", "data", "summary", "gene_MDD_risk_cell_fraction_interaction.csv"))

## Summarise number of significant
(interaction_summary_mdd <- tensorqtl_adj_mdd %>%
    summarise(FDR_ct_sig = sum(FDR_gi < 0.05)) %>%
    left_join(tensorqtl_adj_mdd %>% summarise(bonf_ct_sig = sum(p.bonf_gi < 0.05))))

write_csv(interaction_summary_mdd, file = here("eqtl", "data", "summary", "gene_MDD_risk_cell_fraction_interaction_summary.csv"))

#### Prep MDD Data for Plotting ####
## Merge data
tensorqtl_adj_mdd_anno <- tensorqtl_adj_mdd %>%
    filter(FDR_gi < 0.05) %>%
    dplyr::group_by(BrainRegion) %>%
    arrange(FDR_gi) %>%
    mutate(
        eqtl = paste0(row_number(), ". ", gencodeID, " - ", variant_id),
        eqtl_anno = paste("Gene:", Symbol, "\nSNP:", RS, "\nFDR= ", scales::scientific(FDR_gi, didgits = 3))
    )

tensorqtl_adj_mdd_anno$eqtl <- factor(tensorqtl_adj_mdd_anno$eqtl)
tensorqtl_adj_mdd_anno$eqtl <- fct_reorder(tensorqtl_adj_mdd_anno$eqtl, tensorqtl_adj_mdd_anno$FDR_gi, min)

head(levels(tensorqtl_adj_mdd_anno$eqtl))
tensorqtl_adj_mdd_anno %>% count(eqtl)

express_geno_mdd_cf <- tensorqtl_adj_mdd_anno %>%
    left_join(geno_long_mdd, by = "variant_id") %>%
    left_join(cell_fractions, by = c("BrainRegion", "cell_type", "genoSample")) %>%
    left_join(gene_resid_long, by = c("gencodeID", "genoSample"))

express_geno_mdd_anno <- tensorqtl_adj_mdd_anno %>%
    select(BrainRegion, variant_id, eqtl, eqtl_anno, cell_type) %>%
    left_join(risk_snps_mdd_detials) %>%
    mutate(eqtl_anno = paste0(eqtl_anno, "\n", risk_anno))

#### PLOT ####
## if n/4 %/%
# 2 sACC  61 * will create a page with 1 plot ERROR!
tensorqtl_adj_mdd_anno %>% 
    count(BrainRegion) %>%
    mutate(one_plot = n %% 4 == 1)

# BrainRegion     n one_plot
# <chr>          <int> <lgl>   
# 1 Amygdala       52 FALSE   
# 2 sACC          118 FALSE 

walk(c("Amygdala", "sACC"), function(r) {
    message(r)
    express_cf_sactter <- express_geno_mdd_cf %>%
        filter(BrainRegion == r) %>%
        ggplot(aes(x = cell_fraction, y = expression)) +
        geom_point(aes(color = Genotype), alpha = 0.5, size = .5) +
        geom_smooth(aes(color = Genotype, fill = Genotype), method = lm)

    required_n_pages <- n_pages(express_cf_sactter + facet_wrap_paginate(~ eqtl + cell_type, scales = "free", nrow = 2, ncol = 2))
    message(required_n_pages)

    pdf(here("eqtl", "plots", paste0("eqtl_risk_MDD_cell_fraction_", r, ".pdf")), width = 10)
    for (i in 1:required_n_pages) {
        print(express_cf_sactter +
            facet_wrap_paginate(~ eqtl + cell_type, scales = "free", nrow = 2, ncol = 2, page = i)
            +
            geom_text(
                data = filter(express_geno_mdd_anno, BrainRegion == r),
                aes(label = eqtl_anno), x = -Inf, y = Inf, vjust = "inward", hjust = "inward", size = 3
            ))
    }
    dev.off()
})
