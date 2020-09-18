## read packages
library(jaffelab)
library(readr)
library(SummarizedExperiment)
library(stringr)
library(GenomicRanges)
library(sessioninfo)

## load data
load("../exprs_cutoff/rse_gene.Rdata")
pd = colData(rse_gene)

############# 
# fam file ##

### read in fam
bfile="topmed_mdd_602sample_090120_maf005"
fam = read.table(paste0(bfile, ".fam"), as.is=TRUE)
fam_samples <- paste0(fam$V1, "_", fam$V2)

message("All ",length(unique(pd$genoSample)) , " samples present: ", all(unique(pd$genoSample) %in% fam_samples))

famOut = paste0(fam$V1, " ", fam$V2)[fam_samples %in% pd$genoSample]
cat(famOut, file = "samples_to_extract.txt", sep = "\n")

#### overall extraction
newbfile = "goesHyde_mdd_Genotypes_maf01_geno10_hwe1e6"

## extract
message("\n***** Plink Extract ***** ", Sys.time())
system(paste("plink --bfile", bfile,
             "--keep samples_to_extract.txt --geno 0.1 --maf 0.01 --hwe 0.000001 --make-bed --out", newbfile))

# ## independent and cluster
message("\n***** Plink independent and cluster *****")
system(paste("plink --bfile", newbfile, "--maf 0.1 --indep 100 10 1.25 --out", newbfile))

## MDS components
message("\n***** Plink MDS components *****")
system(paste0("plink --bfile ", newbfile,
              " --cluster --mds-plot 10 --extract ",newbfile, ".prune.in --out ", newbfile))

# ## A transpose
message("\n***** Plink A transpose *****")
system(paste("plink --bfile", newbfile,
             "--recode A-transpose --out", newbfile))

message("\n***** Done Plink ***** ", Sys.time())

# sgejobs::job_single('pull_genotype_data_plink', create_shell = TRUE, queue= 'bluejay', memory = '150G', command = "Rscript pull_genotype_data_plink.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
