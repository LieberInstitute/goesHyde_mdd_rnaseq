library("tidyverse")
library("jaffelab")
library("miniparquet")
library("sessioninfo")
library("here")


## loop over regions
regions <- c(amyg = "Amygdala", sacc = "sACC")

parquet_files <- map(regions,
  ~ list.files(
    path = here("eqtl", "data", "tensorQTL_out", "risk_nominal_mdd"),
    pattern = paste0("*",.x,".*\\.parquet$"),
    full.names = TRUE)
)

map_int(parquet_files, length)
# amyg sacc 
# 230  230

## loci 
loci_list <- read.csv(here("eqtl","data","risk_snps","MDD_GWSig_loci_location_hg38.csv"))

## test data ##
# test_parquet <- parquet_read("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_out/risk_nominal_mdd/gene_sACC_Tcell.cis_qtl_pairs.chr1.parquet")%>%
#   mutate(FDR = p.adjust(pval_gi, "fdr"),
#          chr = ss(variant_id,":"),
#          bp = as.integer(ss(variant_id,":", 2))) 
# 
# 
# head(test_parquet)
# 
# loci_test <- loci_list[1:20,]
# loci_test <- data.frame(locus = 1, chr = "chr1", start = 54489,end = 54500)
# nrow(test_parquet)
# # [1] 4414001
# 
# test_parquet %>% count(chr)

# filter_loci <- function(eqtl_out, loci_tab){
#   
#   eqtl_loci <- pmap(loci_tab, function(locus, chr, start,end){
#     
#     c <- chr ## asignment bug
#     eqtl_out <- eqtl_out %>% filter(chr == c & bp > start, bp < end) %>%
#       mutate(locus = locus)
#     
#     message("locus:",locus," ", chr, " ",start,":",end, " - n = ",nrow(eqtl_out))
#     return(eqtl_out)
#   })
#   
#   eqtl_loci <- do.call("rbind", eqtl_loci)
#   return(eqtl_loci)
# }

# filter_test <- filter_loci(test_parquet, loci_list)
# filter_test  %>%
#   count(locus)
# nrow(filter_test)

read_adj_filter_loci <- function(parquet_files, loci_tab) {
    message("reading ", length(parquet_files)," files - ", Sys.time())
    eqtl_out <- do.call("rbind", map(parquet_files, parquet_read)) %>%
        mutate(FDR = p.adjust(pval_gi, "fdr"),
               chr = ss(variant_id,":"),
               bp = as.integer(ss(variant_id,":", 2))) 
    
    message("unfiltered n pairs: ", nrow(eqtl_out))
    
    ## filter 
    ## loop over loci & label too
    eqtl_loci <- pmap(loci_tab, function(locus, chr, start,end){
      
      c <- chr ## asignment bug
      eqtl_out <- eqtl_out %>% filter(chr == c & bp > start, bp < end) %>%
        mutate(locus = locus)
      
      message("locus:",locus," ", chr, " ",start,":",end, " n = ",nrow(eqtl_out))
      return(eqtl_out)
    })
    
    eqtl_loci <- do.call("rbind", eqtl_loci)
    message("filtered n pairs: ", nrow(eqtl_loci))
    return(eqtl_loci)
}

## test
# eqtl_out_test <- read_adj_filter_loci(parquet_files$amyg$mdd[1:2], loci_list)
# table(eqtl_out_test$chr)


## use pval_gi is p-value for interaction
read_adj_filter_ct <- function(parquet_files, start, end, chr){
    ## get ct from file name
    cell_type <- ss(ss(basename(parquet_files),"\\."),"_",3)
    cell_type_i <- splitit(cell_type)
    ## loop read + get FDR by cell type
    eqtl_out <- map2(cell_type_i, names(cell_type_i),function(cti, ctn){
        message(ctn)
        cell_type_parquets <- parquet_files[cti]
        ct_eqtl_out <- read_adj_filter_loci(cell_type_parquets, loci_tab = loci_list)
        return(ct_eqtl_out %>% mutate(cell_type = ctn, .before = "phenotype_id"))
    })
    ## combine into one table
    eqtl_out_all_ct <- do.call("rbind",eqtl_out)
    rownames(eqtl_out_all_ct) <- NULL
    return(eqtl_out_all_ct)
}


eqtl_out_filtered <- map2(parquet_files, names(parquet_files), function(parq_region, name_region) {
        message(paste("Reading:", name_region, "- "), Sys.time())

        eqtl_out_filtered <- read_adj_filter_ct(parquet_files = parq_region)
        
        fn = here("eqtl", "data", "tensorQTL_loci", paste0("interaction_gene_",name_region, "_loci102.csv"))
        
        message("Saving: ", fn)

        write.csv(eqtl_out_filtered, file = fn, row.names = FALSE)
})

# sgejobs::job_single('filter_interaction_locus', create_shell = TRUE, memory = '50G', command = "Rscript filter_interaction_locus.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
