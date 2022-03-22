library(SummarizedExperiment)
library(tidyverse)
library(jaffelab)
library(sessioninfo)
library(here)

# source(here("eqtl", "code", "utils.R"))
# load(here("data", "MDD_colors.Rdata"), verbose = TRUE)

#### load expression data ####
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)

rd <- as.data.frame(rowData(rse_gene)) %>% select(gencodeID, Symbol)
regions <- c(amyg = "Amygdala", sacc = "sACC")

#### genomewide results ####
## load tables
message("Loading tensorQTL cis results")

eqtl_out <- map(regions,
                ~ read.csv(here("eqtl", "data", "tensorQTL_out", "genomewide_independent", paste0("independent_gene_", .x, ".csv")),
                           row.names = 1) %>%
                    mutate(FDR = p.adjust(pval_perm, "fdr")) # check about this p val -> pval_perm?
)


head(eqtl_out$amyg)
summary(eqtl_out$amyg$FDR)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0003007 0.0003007 0.0011720 0.0051624 0.0053224 0.0412959 
summary(eqtl_out$amyg$pval_perm)

map_int(eqtl_out, nrow)
# amyg sacc 
# 2225 2314

## Not just one results per gene
map_int(eqtl_out, ~ length(unique(.x$phenotype_id)))
# amyg  sacc
# 25085 25085

map_int(eqtl_out, ~ length(unique(.x$variant_id)))
# amyg sacc 
# 1683 1773

#### Extract and save FDR < 0.01 #### 
eqtl_FDR01 <- map(eqtl_out, ~.x %>% filter(FDR < 0.01))
map_int(eqtl_FDR01, nrow)  ## All Signif?
# amyg sacc 
# 1805 1880 

walk2(eqtl_FDR01, names(eqtl_FDR01), 
      ~write.csv(.x, here("eqtl", "data", "tensorQTL_FDR01", "genomewide_independent", paste0("independent_gene_", .y, "_FDR01.csv"))))

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

## Build one table for the data
eqtl_ind <- do.call("rbind", eqtl_out) %>%
    rename(gencodeID = phenotype_id) %>%
    left_join(rd) %>%
    left_join(mdd_snps2) %>%
    replace_na(list(MDD_riskSNP = FALSE))

head(eqtl_ind)

## How many significant?
eqtl_ind %>%
    filter(FDR < 0.01) %>%
    count(BrainRegion)
#   BrainRegion    n
# 1    Amygdala 1805
# 2        sACC 1880

#### How many are risk SNPs? ####
eqtl_ind %>%
    count(MDD_riskSNP)
#   MDD_riskSNP    n
# 1       FALSE 4532
# 2        TRUE    7

eqtl_ind %>%
    filter(FDR < 0.01, MDD_riskSNP)


## Build summary
(ind_summary <- eqtl_ind %>%
        group_by(BrainRegion) %>%
        summarize(n_tested = n(),
                  n_FDR01 = sum(FDR < 0.01),
                  risk_SNPs_tested = sum(MDD_riskSNP),
                  risk_SNPs_FDR01 = sum(MDD_riskSNP & FDR < 0.01))
)
#   BrainRegion n_tested n_FDR01 risk_SNPs_tested risk_SNPs_FDR01
# 1 Amygdala        2225    1805                3               2
# 2 sACC            2314    1880                4               4

write.csv(ind_summary, file = here("eqtl", "data", "summary", "genomewide_independent_summary.csv"))
