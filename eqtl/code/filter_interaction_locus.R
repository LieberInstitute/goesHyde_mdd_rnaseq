library("tidyverse")
library("jaffelab")
library("miniparquet")
library("sessioninfo")
library("here")


## loop over dx
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

test_parquet <- parquet_read("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_out/risk_nominal_mdd/gene_sACC_Tcell.cis_qtl_pairs.chr1.parquet")
test_parquet %>%
  mutate(chr = ss(variant_id,":"),
         bp = as.integer(ss(variant_id,":", 2))) %>%
  tail()


read_adj_filter <- function(parquet_files, start, end, c) {
    message("reading ", length(parquet_files)," files - ", Sys.time())
    eqtl_out <- do.call("rbind", map(parquet_files, parquet_read)) %>%
        mutate(FDR = p.adjust(pval_gi, "fdr"),
               chr = ss(variant_id,":"),
               bp = as.integer(ss(variant_id,":", 2))) 
    
    message("n pairs: ", nrow(eqtl_out))
    
    ## filter 
    ## loop over loci & label too
    eqtl_out <- eqtl_out %>%
        filter(bp > start & bp < end & chr == c)
    
    message("n pairs in locus: ", nrow(eqtl_out))
    # significant_snps <- c(significant_snps, eqt_outl$variant_id)
    return(eqtl_out)
}

## test
eqtl_out_test <- read_adj_filter(parquet_files$amyg$mdd[1:2], start = 7929242, end = 8929242, c = "chr1")
summary(eqtl_out_test$bp)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8000381 8220775 8443875 8486633 8780243 8999779
table(eqtl_out_test$chr)


## loci 
loci_list <- read.csv(here("eqtl","data","risk_snps","MDD_GWSig_loci_location_hg38.csv"))

## use pval_gi is p-value for interaction
read_adj_filter_ct <- function(parquet_files, start, end, chr){
    ## get ct from file name
    cell_type <- ss(ss(basename(parquet_files),"\\."),"_",3)
    cell_type_i <- splitit(cell_type)
    ## loop read + get FDR by cell type
    eqtl_out <- map2(cell_type_i, names(cell_type_i),function(cti, ctn){
        message(ctn)
        cell_type_parquets <- parquet_files[cti]
        ct_eqtl_out <- read_adj_filter(cell_type_parquets, start, end, c = chr)
        return(ct_eqtl_out %>% mutate(cell_type = ctn, .before = "phenotype_id"))
    })
    ## combine into one table
    eqtl_out_all_ct <- do.call("rbind",eqtl_out)
    rownames(eqtl_out_all_ct) <- NULL
    return(eqtl_out_all_ct)
}


 eqtl_out_test_ct <- read_adj_filter_ct(parquet_files$amyg$mdd, start = 8e6, end = 8.5e6, chr = "chr1")

 # head(eqtl_out_test)

eqtl_out_test %>%
    group_by(cell_type) %>%
    summarize(n_pairs = n(),
              n_FDR05 = sum(FDR < 0.05))

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
