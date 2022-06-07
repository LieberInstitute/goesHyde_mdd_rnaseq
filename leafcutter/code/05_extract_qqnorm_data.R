library("SummarizedExperiment")
library("here")
library("jaffelab")
library("purrr")
library("sessioninfo")

load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)
chr_fn <- list.files(here("leafcutter", "data", "clusters"), pattern = "qqnorm_chr.*\\.gz$",
                     full.names = TRUE) 

qqnorm_list <- map(chr_fn, ~read.table(.x, header=T,comment.char="?"))
map_int(qqnorm_list, ncol)
map_int(qqnorm_list, nrow)

## combine in to one table
qqnorm_all <- do.call("rbind", qqnorm_list)
dim(qqnorm_all)
# [1] 210544   1095

corner(qqnorm_all)

colnames(qqnorm_all)[[1]] <- "Chr" 

colnames(qqnorm_all) <- ss(colnames(qqnorm_all), "_")
corner(qqnorm_all)

## fix positions
pos <- qqnorm_all[,c("Chr", "start", "end", "ID")]
# fix end = start + 1
pos$end <- pos$start + 1
head(pos)

## Use genoSample as col names to match VCF
counts <- qqnorm_all[,rse_gene$RNum]
dim(counts)
all(rse_gene$RNum == colnames(counts))
colnames(counts) <- rse_gene$genoSample

corner(counts)

## separate regions and save
region_idx <- splitit(rse_gene$BrainRegion)
names(region_idx)

walk2(region_idx, names(region_idx), function(idx, name){
  
  region_counts <- counts[,idx]
  region_qqnorm <- cbind(pos, region_counts)

  fn <- here("leafcutter", "data", "qqnorm", paste0("qqnorm_",name,".bed"))
  message("saving ", fn)
  
  write.table(region_qqnorm, file= fn,
              quote=F, sep="\t", col.names=T, row.names=F)
  message("gzip...")
  system(paste("gzip", fn))
  
  
})

# sgejobs::job_single('05_extract_qqnorm_data', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript 05_extract_qqnorm_data.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

