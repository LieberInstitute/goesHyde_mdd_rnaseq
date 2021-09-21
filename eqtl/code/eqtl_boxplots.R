
library(tidyverse)
library(VariantAnnotation)
library(jaffelab)
library(here)
library(sessioninfo)

get_resid_expres <- function(rse, mod, tpm = FALSE){
  mod <- mod[colnames(rse),]
  if(tpm) rpkm <- assays(rse)$tpm else rpkm <- recount::getRPKM(rse, "Length")
  ## residualize expression
  exprs <- log2(rpkm + 1)
  exprs <- cleaningY(exprs, mod, P = 1)
  colnames(exprs) <- rse$genoSample
  return(exprs)
}

## load
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)


## define data types
regions <- c(amyg = "Amygdala", sacc = "sACC")
features <- c("gene", "exon", "jxn", "tx")
names(features)  <- features
rse_gene_split <- map(regions, ~rse_gene[,rse_gene$BrainRegion == .x])

## build model
pd <- colData(rse_gene)
mod <- model.matrix(~PrimaryDx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, data = pd)
colnames(mod)

## get residulized expression
resid_expres <- map(rse_gene_split, ~get_resid_expres(.x, mod))

## load eqtl results
eqtl_results <- map(features, function(feat){
  map(regions, function(region) {
    fn = paste0(here("eqtl", "data", "tensorQTL_out",feat),"_", region, "_cis_risk_bpd.csv")
    if(file.exists(fn)){
      csv = read_csv(fn, show_col_types = FALSE)
      return(csv)
    } else{
      return(NULL)
    }
  })
})

map(eqtl_results, ~map(.x, dim))
eqtl_results$gene$amyg %>% dplyr::count(pval_beta < 0.05)

gene_amyg <- read_csv(output_csv[[1]], show_col_types = FALSE)

## get genotypes
vcf_risk_bpd <- readVcf(here("eqtl","data","risk_snps", "LIBD_maf01_gwas_BPD.vcf.gz"))

vcf_split <- map(rse_gene_split, ~vcf_risk_bpd[,.x$genoSample])
map(vcf_split, dim)

geno_long <- map(vcf_split, ~as.data.frame(geno(.x)$GT) %>% 
                   rownames_to_column("variant_id") %>%
                   pivot_longer(!variant_id, names_to = "genoSample", values_to = "Genotype") %>%
                   mutate(Genotype = case_when(Genotype == "0|0"~"0",
                                               Genotype == "1|1"~"2",
                                               TRUE~"1")))

## Join eqtl results, residual espression, and genotype

eqtl_filter <- eqtl_results$gene$amyg %>% arrange(pval_beta)%>% dplyr::select(phenotype_id, variant_id, pval_beta) %>% dplyr::slice(1:10)

eqtl_filter <- map(eqtl_results["gene"], ~map(.x, function(e){
  e <- e %>% 
    arrange(pval_beta) %>% 
    dplyr::select(phenotype_id, variant_id, pval_beta) %>% 
    dplyr::slice(1:10) %>%
    mutate(eqtl = paste0("Gene: ", phenotype_id, "\n", "SNP: ", variant_id))
  
  e$eqtl <- factor(e$eqtl)
  e$eqtl <- fct_reorder(e$eqtl, e$pval_beta, min)
  return(e)
  } 
))

levels(eqtl_filter$gene$amyg$eqtl)

eqtl_expres_geno <- pmap(list(ef = eqtl_filter, gl = geno_long, re = resid_expres),
                         function(ef, gl, re){
                           map(ef, ~re[.x$phenotype_id,] %>%
                               as.data.frame() %>%
                               rownames_to_column("phenotype_id") %>%
                               pivot_longer(!phenotype_id, names_to = "genoSample", values_to = "Expression") %>%
                               left_join(.x %>% dplyr::select(phenotype_id, variant_id, pval_beta, eqtl), by = "phenotype_id") %>%
                               left_join(gl, by = c("genoSample", "variant_id"))
                               )
                         })

walk2(eqtl_expres_geno, names(eqtl_expres_geno), function(feat_data, feat_name){
  walk2(feat_data, feat_name, function(region_data, region_name){
    message(here("eqtl","plots",paste0("eqtl_boxplot_bpd_risk_", feat_name, "_", region_name,".png")))
    eqtl_boxplot <- ggplot(region_data, aes(Genotype, Expression, fill = Genotype)) +
      geom_boxplot() +
      facet_wrap(~eqtl, nrow = 2) + 
      geom_text(aes(x = 1, y = 4.75, label = paste("p =", format.pval(pval_beta))),
                vjust = "inward", hjust = "inward", size = 3) +
      labs(y = "Residualized Expression", title = "cis-eQTLs BPD risk SNPs",
           subtitle = paste(feat_name, region_name)) +
      guides(fill="none") +
      theme_bw()
    
    ggsave(eqtl_boxplot, filename = here("eqtl","plots",paste0("eqtl_boxplot_bpd_risk_", feat_name, "_", region_name,".png")), width = 12)
    
  })
} )


