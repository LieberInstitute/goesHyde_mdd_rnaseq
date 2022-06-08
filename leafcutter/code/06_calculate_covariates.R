library("SummarizedExperiment")
library("jaffelab")
library("purrr")
library("sva")
library("here")
library("sessioninfo")

## get region
args <- commandArgs(trailingOnly = TRUE)
region <- args[[1]]
message("region = ",region)
#### load data ####
load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)

rse_gene <- rse_gene[,rse_gene$BrainRegion == region]
pd = as.data.frame(colData(rse_gene))


covar_format <- function(data, rn) {
  data <- as.data.frame(data)
  rownames(data) <- rn
  data <- t(data)
  data <- as.data.frame(data) %>% tibble::rownames_to_column("id")
  return(data)
}

## define model
mod <- model.matrix(~ PrimaryDx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, data = pd)

fn <- here("leafcutter", "data", "qqnorm", paste0("qqnorm_",region,".bed.gz"))

message("reading ", basename(fn), " - ", Sys.time())
splice <- read.table(fn,sep='\t',comment.char='$',header=T)

## drop meta data columns
splice <- splice[,5:ncol(splice)]
  
## run PCA
message("running PCA - ", Sys.time())
sPCA <- prcomp(t(splice))

message("running num.sv - ", Sys.time())
numPCs <- num.sv(splice, mod)

message("Done - ", Sys.time())
message("k = ", numPCs)
sPCs <- sPCA$x[,1:numPCs]


covsSplice <- cbind(mod[,-1],sPCs)
covsSplice <- covar_format(covsSplice, rse_gene$genoSample)

write.table(covsSplice,file= here("leafcutter", "data","covariates",paste0("covariates_",region,".txt")),
            sep="\t",quote=F,col.names=T,row.names=F)

# sgejobs::job_single('06_calculate_covariates', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript 06_calculate_covariates.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
