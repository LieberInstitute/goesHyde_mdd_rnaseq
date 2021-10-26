library(tidyverse)
library(jaffelab)
library(miniparquet)
library(sessioninfo)
library(here)

features <- c("gene", "exon", "jxn", "tx")
names(features) <- features

regions <- c(amyg = "Amygdala", sacc = "sACC")

parquet_files <- map(features,
                     function(f) map(regions,
                                     ~list.files(path = here("eqtl", "data", "tensorQTL_out","cis_genomewide_nominal"),
                                                 pattern = paste0(f,'_',.x),
                                                 full.names = TRUE)))

read_adj_filter <- function(parquet_files, cutoff = 0.01){
  
  eqtl_out <- do.call("rbind", map(parquet_files, parquet_read)) %>%
    mutate(FDR = p.adjust(pval_nominal, "fdr")) %>%
    filter(FDR < cutoff)
  
  return(eqtl_out)
}


eqtl_out_filtered <- map2(parquet_files, names(parquet_files), function(parq_feat, names_feat){

  map2(parq_feat, names(parq_feat), function(parq_region, names_region){
    
    message(paste("Reading:",names_feat, names_region, "- "),Sys.time())
    
    # eqtl_out_filtered <- read_adj_filter(parquet_files = parq_region[1])
    return(here("eqtl", "data", "tensorQTL_out","cis_genomewide_nominal",
         paste0(names_feat, "_",names_region,"_FDR01.csv")))
    # write_csv(eqtl_out_filtered, file = here("eqtl", "data", "tensorQTL_out","cis_genomewide_nominal",
    #                                          paste0(names_feat, "_",names_region,"_FDR01.csv")))
  })
})

#sgejobs::job_single('filter_genomewide_eqtl', create_shell = TRUE, memory = '100G', command = "Rscript filter_genomewide_eqtl.R")
## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
