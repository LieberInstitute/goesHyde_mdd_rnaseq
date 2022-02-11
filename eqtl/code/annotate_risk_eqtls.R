library(SummarizedExperiment)
library(tidyverse)
library(here)

## Load rse objects
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_exon.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_jxn.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_tx.Rdata"), verbose = TRUE)

## get rowData
rd_cols <- c("gencodeID", "Symbol", "gene_type", "Class")

rd <- map(c(gene = rse_gene, exon = rse_exon, jxn = rse_jxn, tx = rse_tx),
          ~rowData(.x)  %>%
            as.data.frame() %>%
            rownames_to_column("phenotype_id")
)

gene_details <- rd$gene %>% 
  select(gencodeID, Symbol, gene_type, Class)

## Fix some rownames
rd$gene <- rd$gene %>% 
  select(phenotype_id, gencodeID, Symbol, gene_type, Class)
  
rd$exon <- rd$exon %>% 
  select(phenotype_id, gencodeID, Symbol, gene_type, Class)

rd$jxn <- rd$jxn %>%
  select(phenotype_id, gencodeID = gencodeGeneID) %>%
  left_join(gene_details)

rd$tx <- rd$tx %>%
  select(phenotype_id, gencodeID = gene_id) %>%
  left_join(gene_details)

map(rd, head)

risk_eqtl_anno <- left_join(risk_eqtl, rd)
head(risk_eqtl_anno)

risk_eqtl_anno %>% 
  filter(FDR < 0.05) %>% 
  count(Symbol) %>%
  arrange(-n)

## Risk snp data
# vcf_fn <- map(c(MDD = 'MDD', BPD = "BPD"),
#               ~here("eqtl", "data", "risk_snps", paste0("LIBD_maf01_gwas_",.x,".vcf.gz")))
# 
# bpd <- VariantAnnotation::readVcf(vcf_fn$BPD)
# 
# snp_info <- map(vcf_fn, ~VariantAnnotation::readInfo(.x, x = "DP"))
# map(snp_info, head)

## Risk eQTL files

walk2(rd, names(rd), function(feat_rd, name){
  fn <- list.files(path = here("eqtl", "data", "risk_snps_eqtl"), 
                   full.names = TRUE, 
                   pattern = paste0("*", name,"*"))
  
  walk(fn, function(f){
    message(fn)
    e <- read.csv(f) %>%
      left_join(feat_rd)
    write.csv(e, file = f)
    
  })
})

