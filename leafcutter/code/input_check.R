
library("SummarizedExperiment")
library("here")

load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)

## Get all 1091 the junction count file paths to convert to .junc files
junc_paths <- gsub("accepted_hits.sorted.bam","junctions_primaryOnly_regtools.count", gsub("HISAT2_out","Counts/junction", rse_gene$bamFile))
head(junc_paths)

length(junc_paths)
# [1] 1091

all(file.exists(junc_paths))
# [1] TRUE


cat(junc_paths, file = here("leafcutter","data","junc_count_files.txt"), sep = "\n")
