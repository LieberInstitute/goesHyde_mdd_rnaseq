library("SummarizedExperiment")
library("tidyverse")
library("jaffelab")
library("sessioninfo")
library("here")

#### Check nums from log files ####
nominal_log_fn <- list.files(here("leafcutter","code", "logs"), pattern = "07_run_tensorqtl", full.names = TRUE)

nominal_log <- map(nominal_log_fn, readLines)
pair <- map_chr(nominal_log, ~.x[[grep("Reading Expression files:",.x)]])
names(nominal_log) <- ss(ss(basename(pair),"_",2),"\\.")

time_lines <- map(nominal_log, ~.x[grep("time elapsed: ",.x)])
phenotype <- map_dbl(nominal_log, ~parse_number(.x[grep("\\^* \\d+ phenotypes",.x)][[1]]))
runtime <- map_dbl(nominal_log, ~max(parse_number(.x[grep("time elapsed: ",.x)])))

filter_log <- readLines(here("leafcutter","code", "logs", "08_filter_sQTL.txt"))

n_pairs <- parse_number(filter_log[grep("n pairs:", filter_log)])
n_pairs_FDR05 <- parse_number(gsub("n pairs FDR<0.05:","",filter_log[grep("n pairs FDR<0.05:", filter_log)]))

nominal_log_data <- tibble(feat ="Splice",
                           region = names(runtime),
                           n_feat = phenotype,
                           runtime,
                           n_pairs,
                           n_pairs_FDR05 )
# 
# # A tibble: 2 × 5
#   data     n_feat runtime   n_pairs n_pairs_FDR05
# <chr>     <dbl>   <dbl>     <dbl>         <dbl>
# 1 Amygdala 209476    19.2 451404371       7851177
# 2 sACC     209476    19.2 449487451       8324138

write.csv(nominal_log_data, file = here("leafcutter","data", "LC_sQTL_summary.csv"))

#### Check out tables ####
sqtl_fn <- list.files(here("leafcutter","data","tensorQTL_FDR05"), full.names = TRUE)

sQTL <- read.csv(sqtl_fn[[1]])

length(unique(sQTL$phenotype_id))
# [1] 91292
length(unique(sQTL$variant_id))
# [1] 1633844

pheno_count <- sQTL %>% group_by(phenotype_id) %>% count()

summary(pheno_count$n)

pheno_count %>% arrange(-n)

#### Check phenotype ####
load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)


cluster_roi <- function(clust_name){
  clust_name <- unlist(strsplit(clust_name, ":"))
  roi <- GRanges(seqnames=clust_name[[1]], ranges=as.integer(clust_name[[2]]):as.integer(clust_name[[3]]))
  return(roi)
}

cluster_roi("chr1:150269330:150270518:clu_17903_NA")

subsetByOverlaps(rse_gene, cluster_roi("chr1:150269330:150270518:clu_17903_NA"))

roi <- GRanges(seqnames="chr1", ranges=600505:630761)
rowData(subsetByOverlaps(rse_gene, roi))


# DataFrame with 2 rows and 10 columns
# Length         gencodeID       ensemblID              gene_type      Symbol  EntrezID       Class
# <integer>       <character>     <character>            <character> <character> <integer> <character>
#   ENSG00000225972.1       372 ENSG00000225972.1 ENSG00000225972 unprocessed_pseudogene    MTND1P23        NA       InGen
# ENSG00000225630.1      1044 ENSG00000225630.1 ENSG00000225630 unprocessed_pseudogene    MTND2P28        NA       InGen
# NumTx meanExprs passExprsCut
# <integer> <numeric>    <logical>
#   ENSG00000225972.1         1    5.4378         TRUE
# ENSG00000225630.1         1   46.0940         TRUE
