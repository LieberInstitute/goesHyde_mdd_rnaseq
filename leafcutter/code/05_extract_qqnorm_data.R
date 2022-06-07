library("SummarizedExperiment")
library("here")
library("jaffelab")

load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)

# all <- read.table(here("leafcutter", "data", "clusters", "leafcutter_perind.counts.gz"),
#                   header=T,comment.char="?", row.names= 1)

chr_fn <- list.files(here("leafcutter", "data", "clusters"), pattern = "qqnorm_chr*.gz") 

all <- read.table(here("leafcutter", "data", "clusters", "leafcutter_perind.counts.gz.qqnorm_chr1.gz"),
                  header=T,comment.char="?")

colnames(all)[[1]] <- "Chr" 
dim(all)
# [1] 289778   1091

colnames(all) <- ss(colnames(all), "_")
corner(all)

## fix end = start + 1
all$end <- all$start + 1

## Use genoSample as col names to match VCF
pos <- all[,c("Chr", "start", "end", "ID")]
head(pos)

counts <- all[,rse_gene$RNum]
all(rse_gene$RNum == colnames(counts))
colnames(counts) <- rse_gene$genoSample

corner(counts)

## separate regions and save
region_idx <- splitit(rse_gene$BrainRegion)
names(region_idx)

amy <- counts[,region_idx$Amygdala]
amy <- cbind(pos, counts)
corner(amy)
# [1]  10 540
amy_fn <- here("leafcutter", "data", "qqnorm", "qqnorm_Amygdala_chr1.bed")
write.table(amy, file= amy_fn,
            quote=F, sep="\t", col.names=T, row.names=F)
system(paste("gzip", amy_fn))

sacc <- all[,region_idx$sACC]
dim(sacc)
# [1]  10 551
sacc_fn <- here("leafcutter", "data", "qqnorm", "qqnorm_sACC.txt")
write.table(sacc, file = sacc_fn,
            quote=F, sep="\t", col.names=T, row.names=F)
system(paste("gzip", sacc_fn))

# sgejobs::job_single('extract_qqnorm_data', create_shell = TRUE, queue= 'bluejay', memory = '5G', command = "Rscript extract_qqnorm_data.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

