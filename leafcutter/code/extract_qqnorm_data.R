library("SummarizedExperiment")
library("here")
library("jaffelab")

load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)

all <- read.table(here("leafcutter", "data", "clusters", "leafcutter_perind.counts.gz.PCs"),
                  header=T,comment.char="?")

colnames(all) <- ss(colnames(all), "_")

## Use genoSample as col names to match VCF
all <- all[,rse_gene$RNum]
all(rse_gene$RNum == colnames(all))
colnames(all) <- rse_gene$genoSample

region_idx <- splitit(rse_gene$BrainRegion)
names(region_idx)

## separate regions and save
amy <- all[,br_amy]
write.table(amy, file="/users/schadinh/lieber/qqnorm/amygdala_qqnorm_20190329.txt", quote=F, sep="\t", col.names=T, row.names=F)

sacc <- all[,br_sacc]
write.table(sacc, file="/users/schadinh/lieber/qqnorm/sacc_qqnorm_20190329.txt", quote=F, sep="\t", col.names=T, row.names=F)

# sgejobs::job_single('extract_qqnorm_data', create_shell = TRUE, queue= 'bluejay', memory = '5G', command = "Rscript extract_qqnorm_data.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
