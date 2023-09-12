library("tidyverse")
library("jaffelab")
library("miniparquet")
library("sessioninfo")
library("here")

source("utils.R")
 
## loop over dx
# dx <- c("mdd", "bpd")
dx <- c("mdd")
names(dx) <- dx

regions <- c(amyg = "Amygdala", sacc = "sACC")

parquet_files <- map(
    regions,
    function(r) {
        map(
            dx,
            ~ list.files(
                path = here("eqtl", "data", "tensorQTL_out", paste0("risk_nominal_", .x)),
                pattern = paste0("*",r,".*\\.parquet$"),
                full.names = TRUE
            )
        )
    }
)

map_depth(parquet_files, 2, length)

## use pval_gi is p-value for interaction
read_adj_filter_ct <- function(parquet_files, cutoff = 0.05){
    ## get ct from file name
    cell_type <- ss(ss(basename(parquet_files),"\\."),"_",3)
    cell_type_i <- splitit(cell_type)
    ## loop read + get FDR by cell type
    eqtl_out <- map2(cell_type_i, names(cell_type_i),function(cti, ctn){
        message(ctn)
        cell_type_parquets <- parquet_files[cti]
        ct_eqtl_out <- read_adj_filter(cell_type_parquets, cutoff, pval_name = "pval_gi")
        return(ct_eqtl_out %>% mutate(cell_type = ctn, .before = "phenotype_id"))
    })
    ## combine into one table
    eqtl_out_all_ct <- do.call("rbind",eqtl_out)
    rownames(eqtl_out_all_ct) <- NULL
    return(eqtl_out_all_ct)
}

## test
# eqtl_out_test <- read_adj_filter_ct(parquet_files$amyg$mdd)
# head(eqtl_out_test)

# eqtl_out_test %>%
#     group_by(cell_type) %>%
#     summarize(n_pairs = n(),
#               n_FDR05 = sum(FDR < 0.05))

eqtl_out_filtered <- map2(parquet_files, names(parquet_files), function(parq_region, name_region) {
    map2(parq_region, names(parq_region), function(parq_dx, name_dx) {
        message(paste("Reading:", name_region, name_dx, "- "), Sys.time())

        eqtl_out_filtered <- read_adj_filter_ct(parquet_files = parq_dx)
        
        fn = here("eqtl", "data", "tensorQTL_FDR05", paste0("risk_nominal_",name_dx),
          paste0("risk_gene_",name_region, "_",name_dx, "_FDR05.csv"))
        
        message("Saving: ", fn)

        write_csv(eqtl_out_filtered, file = fn)
    })
})
# significant_snps <- unique(significant_snps)
# length(significant_snps)
# cat(significant_snps, sep = "\n", file = here("eqtl", "data", "signif_snps", "significant_snps.txt"))

# sgejobs::job_single('filter_risk_nominal', create_shell = TRUE, memory = '10G', command = "Rscript filter_risk_nominal.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
