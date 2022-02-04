#### Subset nominal results to defined risk SNPs ###

library(tidyverse)
library(here)
library(sessioninfo)
library(miniparquet)

#### Risk SNP Data ####
dx <- c(mdd = "MDD", bpd = "BPD")

risk_SNPs <- map(dx, ~readLines(here("eqtl", "data", "risk_snps", paste0(.x,"_risk_snps.txt"))))
map_int(risk_SNPs, length)
# mdd  bpd 
# 9650   63 
map(risk_SNPs, head)

#### List Parquet Data ####
# features <- c("gene", "exon", "jxn", "tx")
features <- c("exon", "jxn", "tx") ## rerun with more mem w/o gene level
names(features) <- features

regions <- c(amyg = "Amygdala", sacc = "sACC")


parquet_files <- map(
  features,
  function(f) {
    map(
      regions,
      ~ list.files(
        path = here("eqtl", "data", "tensorQTL_out", "genomewide_nominal"),
        pattern = paste0(f, "_", .x),
        full.names = TRUE
      )
    )
  }
)

#### Filter Data ####

walk2(parquet_files, names(parquet_files), function(parq_feat, names_feat) {
  walk2(parq_feat, names(parq_feat), function(parquet_region, names_region) {
    
    message(paste("Reading:", names_feat, names_region, "- "), Sys.time())
    
    ## Read and combine data across CHR
    eqtl_out <- do.call("rbind", map(parquet_region, parquet_read)) %>%
      mutate(FDR = p.adjust(pval_nominal, "fdr"))
    
    message("nrow: ", nrow(eqtl_out))
    
    ## Filter over lists of risk SNPs
    map2(risk_SNPs, dx, function(snps, dx_name){
      message("filtering ", dx_name, " SNPs")
      
      eqtl_out_filtered <- eqtl_out  %>%
        filter(variant_id %in% snps) 
      
      write_csv(eqtl_out_filtered, file = here(
        "eqtl", "data", "risk_snps_eqtl", 
        paste0("nominal_",dx_name,"_riskSNPs_",names_feat, "_", names_region,".csv")
      ))
    })
    
  })
})


# sgejobs::job_single('subset_genomewide_nominal_riskSNPs', create_shell = TRUE, memory = '150G', command = "Rscript subset_genomewide_nominal_riskSNPs.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
