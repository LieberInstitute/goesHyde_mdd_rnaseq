library("tidyverse")
library("jaffelab")
library("miniparquet")
library("sessioninfo")
library("here")

regions <- c(amyg = "Amygdala", sacc = "sACC")
parquet_files <- map(regions, ~list.files(here("leafcutter","data","tensorQTL_out"), pattern = .x, full.names = TRUE))
                     
read_adj_filter <- function(parquet_files, cutoff = 0.05) {
  eqtl_out <- do.call("rbind", map(parquet_files, parquet_read)) %>%
    mutate(FDR = p.adjust(pval_nominal, "fdr")) 
  
  message("n pairs: ", nrow(eqtl_out))
  ## filter 
  eqtl_out <- eqtl_out %>%
    filter(FDR < cutoff)
  message("n pairs FDR<", cutoff, ": ", nrow(eqtl_out))
  return(eqtl_out)
}


eqtl_out_filtered <- map2(parquet_files, regions, function(parq, region) {
    message(paste("Reading:", region, "- "), Sys.time())
    
    eqtl_out_filtered <- read_adj_filter(parquet_files = parq)
    
    file_out <- paste0("LC_nominal_", region, "_FDR05.csv")
    message(paste("writing:", region, "- "), Sys.time())
    
    write_csv(eqtl_out_filtered, file = here("leafcutter", "data", "tensorQTL_FDR05", file_out))
})

# sgejobs::job_single('08_filter_sQTL', create_shell = TRUE, memory = '100G', command = "Rscript 08_filter_sQTL.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
