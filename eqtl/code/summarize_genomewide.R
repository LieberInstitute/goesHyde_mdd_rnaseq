library("SummarizedExperiment")
library("VariantAnnotation")
library("tidyverse")
library("jaffelab")
library("ggforce")
library("ggrepel")
library("sessioninfo")
library("here")

# source(here("eqtl", "code", "utils.R"))
load(here("data", "MDD_colors.Rdata"), verbose = TRUE)

#### load expression data ####
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)

rd <- as.data.frame(rowData(rse_gene)) %>% select(gencodeID, Symbol)
regions <- c(amyg = "Amygdala", sacc = "sACC")

#### Nomninal ####
## check out log files for big numbers 
nominal_log_fn <- list.files(here("eqtl","code", "logs"), pattern = "tensorqtl_genomewide_nominal_gpu*", full.names = TRUE)

nominal_log <- map(nominal_log_fn, readLines)
pair <- map_chr(nominal_log, ~.x[[grep("Reading Expression files:",.x)]])
names(nominal_log) <- ss(ss(pair,"/", 5),"\\.")

time_lines <- map(nominal_log, ~.x[grep("time elapsed: ",.x)])
phenotype <- map_dbl(nominal_log, ~parse_number(.x[grep("\\^* \\d+ phenotypes",.x)][[1]]))
runtime <- map_dbl(nominal_log, ~max(parse_number(.x[grep("time elapsed: ",.x)])))

filter_log <- readLines(here("eqtl","code", "logs", "filter_genomewide_nominal.txt"))

n_pairs <- parse_number(filter_log[grep("n pairs:", filter_log)])
n_pairs_FDR05 <- parse_number(gsub("n pairs FDR<0.05:","",filter_log[grep("n pairs FDR<0.05:", filter_log)]))

nominal_log_data <- tibble(data = names(runtime),
                       n_feat = phenotype,
                       runtime,
                       n_pairs,
                       n_pairs_FDR05 ) %>%
  separate(data,  into = c("feat","region"), remove = FALSE)

# feat   region   n_feat runtime   n_pairs n_pairs_FDR05
# <chr>  <chr>     <dbl>   <dbl>     <dbl>         <dbl>
#   1 gene   Amygdala  25132    2.52  53637238        771221
# 2 gene   sACC      25132    2.55  53403152        853567
# 3 exon   Amygdala 399023   45.6  840508027       5960604
# 4 exon   sACC     399023   45.6  836725998       7192434
# 5 jxn    Amygdala 304125   37.0  633281109       9381590
# 6 jxn    sACC     304125   36.9  630534075      10434226
# 7 tx     Amygdala  79565    8.21 167717620       1684554
# 8 tx     sACC      79565    8.22 166980983       1906150
# 9 Splice Amygdala 209476   19.2  451404371       7851177
# 10 Splice sACC     209476   19.2  449487451       8324138

sQTL_log_data <- read.csv(here("leafcutter","data", "LC_sQTL_summary.csv"), row.names = 1) %>%
  mutate(data = paste0(feat, "_",region))

nominal_log_data <- rbind(nominal_log_data, sQTL_log_data)

nominal_log_data$feat <- factor(nominal_log_data$feat, levels = c("gene", "tx", "exon","jxn","Splice"))

nominal_runtime_scatter <- ggplot(nominal_log_data, aes(x = n_pairs, y = runtime))+
  geom_point() +
  geom_text_repel(aes(label = data)) +
  labs(y = "Runtime (min)") +
  theme_bw()

ggsave(nominal_runtime_scatter, filename = here("eqtl", "plots", "nominal_runtime_scatter.png"))

## checks 
nominal_log_data %>%
  group_by(data) %>%
  transmute(precent_FDR05 = 100*n_pairs_FDR05/n_pairs)
# data          precent_FDR05
# <chr>                 <dbl>
#   1 gene_Amygdala       1.44 
# 2 gene_sACC             1.60 
# 3 exon_Amygdala         0.709
# 4 exon_sACC             0.860
# 5 jxn_Amygdala          1.48 
# 6 jxn_sACC              1.65 
# 7 tx_Amygdala           1.00 
# 8 tx_sACC               1.14


## out of 11M how many SNPs are in pairs for each combo?
nominal_log_data %>%
  group_by(data) %>%
  transmute(mean_snp = n_pairs/n_feat)
# data          mean_snp
# <chr>            <dbl>
# 1 gene_Amygdala    2134.
# 2 gene_sACC        2125.
# 3 exon_Amygdala    2106.
# 4 exon_sACC        2097.
# 5 jxn_Amygdala     2082.
# 6 jxn_sACC         2073.
# 7 tx_Amygdala      2108.
# 8 tx_sACC          2099.

## Significant barplot
n_signif <- nominal_log_data %>%
  mutate(n_pairs_ns = n_pairs - n_pairs_FDR05) %>%
  select(data, n_pairs, n_pairs_ns, n_pairs_FDR05) %>%
  pivot_longer(!c(data, n_pairs), names_to = "Signif", values_to = "n", names_prefix ="n_pairs_")

pairs_barplot <- n_signif %>%
  ggplot(aes(x = data, fill = Signif)) +
  geom_col(aes(y = n)) +
  geom_text(aes(label= ifelse(Signif == "FDR05",
                              formatC(n, format = "e", digits = 2)
                              ,""), 
                y = n_pairs), size = 3)+
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(pairs_barplot, filename = here("eqtl", "plots", "pairs_barplot.png"))

##
signif_barplot <- nominal_log_data  %>%
  ggplot(aes(x = feat, fill = region)) +
  geom_col(aes(y = n_pairs_FDR05), position = "dodge") +
  theme_bw() + 
  scale_fill_manual(values = mdd_BrainRegion_colors) +
  title("Nominal QTL Results")

ggsave(signif_barplot, filename = here("eqtl", "plots", "signif_barplot.png"))


## load nominal files
nom_fn <- list.files(here("eqtl", "data", "tensorQTL_FDR05","genomewide_nominal"), pattern = "*.csv", full.names = TRUE)
names(nom_fn) <- gsub("_FDR05.csv","", basename(nom_fn))

## do the n FDR < 05 match?
# nom_eqtl <- map(nom_fn, read.csv)
# nom_row <- map_int(nom_eqtl, nrow)
# # nominal_exon_amyg nominal_exon_sacc nominal_gene_amyg nominal_gene_sacc  nominal_jxn_amyg  nominal_jxn_sacc 
# # 5599118           5449203            625776            649396           7422175           7935775 
# # nominal_tx_amyg   nominal_tx_sacc 
# # 1365378           1444336 
# 
##yes! 
# nominal_log_data %>%
#   mutate(data_short = gsub("amygdala",'amyg',tolower(data))) %>%
#   left_join(stack(nom_row) %>% rename(data_short = ind)) %>%
#   mutate(n_pairs_FDR05 == values)

nom_eqtl_all <- do.call("rbind", nom_eqtl)
nom_eqtl_all <- nom_eqtl_all %>% 
  rownames_to_column("test") %>%
  separate(test, into = c("test","feature","region"))
  

#### genomewide results ####
## load tables
message("Loading tensorQTL cis results")

eqtl_out <- map(
    regions,
    ~ read.csv(here("eqtl", "data", "tensorQTL_out", "genomewide_cis", paste0("cis_gene_", .x, ".csv"))) %>%
        rename(FDR = qval)
    # %>%
    #     separate(variant_id, into = c("chr", NA, NA, NA), sep = ":", remove = FALSE)
)


head(eqtl_out$amyg)
summary(eqtl_out$amyg$tss_distance)
summary(eqtl_out$amyg$FDR)

map_int(eqtl_out, nrow)
# amyg  sacc
# 25085 25085

## Remember one result per gene = 'phenotype_id'
map_int(eqtl_out, ~ length(unique(.x$phenotype_id)))
# amyg  sacc
# 25085 25085

map_int(eqtl_out, ~ length(unique(.x$variant_id)))
# amyg  sacc
# 16339 16027

#### Extract and save FDR < 0.01 #### 
eqtl_FDR05 <- map(eqtl_out, ~.x %>% filter(FDR < 0.05))
map_int(eqtl_FDR05, nrow)
# amyg sacc 
# 2080 2112

walk2(eqtl_FDR05, names(eqtl_FDR05), 
      ~write.csv(.x, here("eqtl", "data", "tensorQTL_FDR05", "genomewide_cis", paste0("cis_gene_", .y, "_FDR05.csv"))))

#### Extract and save FDR < 0.01 #### 
eqtl_FDR01 <- map(eqtl_out, ~.x %>% filter(FDR < 0.01))
map_int(eqtl_FDR01, nrow)
# amyg sacc 
# 1326 1367

walk2(eqtl_FDR01, names(eqtl_FDR01), 
      ~write.csv(.x, here("eqtl", "data", "tensorQTL_FDR01", "genomewide_cis", paste0("cis_gene_", .y, "_FDR01.csv"))))

#### Combine Region Data and Mutate ####
## combine regions
eqtl_out <- map2(eqtl_out, regions, ~.x %>% mutate(BrainRegion = .y))

## load risk SNP data
mdd_snps <- read.delim(here("eqtl", "data", "risk_snps", "PGC_MDD_genome-wide_significant_Jan2022.txt")) %>%
    filter(!is.na(bp)) %>%
    rename(variant_id = chr_bp_ref_alt)

## Most won't be in the table
mdd_snps2 <- mdd_snps %>%
    select(variant_id) %>%
    mutate(MDD_riskSNP = TRUE)

## Build one table for the cis data
eqtl_cis <- do.call("rbind", eqtl_out) %>%
    rename(gencodeID = phenotype_id) %>%
    left_join(rd) %>%
    left_join(mdd_snps2) %>%
    replace_na(list(MDD_riskSNP = FALSE))

head(eqtl_cis)

## How many significant?
eqtl_cis %>%
    filter(FDR < 0.01) %>%
    count(BrainRegion)
# BrainRegion    n
# 1    Amygdala 1326
# 2        sACC 1367

#### How many are risk SNPs? ####
eqtl_cis %>%
    count(MDD_riskSNP)

eqtl_cis %>%
    filter(FDR < 0.01, MDD_riskSNP)


## Build summary
(cis_summary <- eqtl_cis %>%
    group_by(BrainRegion) %>%
    summarize(total_variants = sum(num_var),
              cis_eqtls = n(),
              cis_FDR05 = sum(FDR < 0.05),
              cis_FDR01 = sum(FDR < 0.01))
)

write.csv(cis_summary, file = here("eqtl", "data", "summary", "genomewide_cis_summary.csv"))

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

map_int(eqtl_out, ~ length(unique(.x$variant_id)))
map_int(eqtl_out, ~ sum(unique(.x$variant_id) %in% rs_id$variant_id))

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
eqtl_out_top$gene$amyg %>%
    count(gencodeID) %>%
    arrange(n) %>%
    left_join(rd)

## annotation for boxplots - order eqtl pairs by FDR
eqtl_out_anno <- map_depth(eqtl_out_top, 2, function(e) {
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
map_depth(eqtl_out_anno, 2, head)
map_depth(eqtl_out_anno, 2, nrow)

## merge genotype and expression data for boxplots
express_geno <- map2(
    eqtl_out_anno, resid_expres_split_long, function(eqtl, expres) {
        map2(eqtl, expres, ~ .x %>%
            inner_join(geno_long, by = "variant_id") %>%
            inner_join(.y, by = c("gencodeID", "genoSample")) %>%
            left_join(genoDx, by = "genoSample"))
    }
)

pwalk(
    list(expres = express_geno, anno = eqtl_out_anno, feat_name = names(express_geno)), function(expres, anno, feat_name) {
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
