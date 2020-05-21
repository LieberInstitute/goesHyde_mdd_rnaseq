library(sessioninfo)
library(here)

# create rdas directory
dir.create(here("eqtl", "genomewide", "rdas"), showWarnings = FALSE)
########################

######################## get overlapping snps
## load SNP data

## Amyg / sACC snps
load(here("genotype_data", "goesHyde_bipolarMdd_Genotypes_n593.rda"), verbose = TRUE) # need ~50G memory to load
snpMap$pos_hg19 <- paste0(snpMap$CHR, ":", snpMap$POS)
snpMap_amyg <- snpMap

## DLPFC snps
load("/dcl01/ajaffe/data/lab/brainseq_phase1/genotype_data/brainseq_phase1_Genotypes_n732.rda", verbose = TRUE) # need 60G ram to load
snpMap$pos_hg19 <- paste0(snpMap$CHR, ":", snpMap$POS)
rm(snp, mds)


ind <- which(snpMap$SNP %in% snpMap_amyg$SNP)
snpMap <- snpMap[ind, ]
snpMap_amyg <- snpMap_amyg[rownames(snpMap), ]

# ## check
stopifnot(identical(snpMap$POS, snpMap_amyg$POS))
stopifnot(identical(snpMap$newRef, snpMap_amyg$newRef))
stopifnot(identical(snpMap$rsNumGuess, snpMap_amyg$rsNumGuess))

snpMapKeep <- snpMap

save(snpMapKeep, file = here("eqtl", "genomewide", "rdas", "overlappingSNPs.rda"))

#sgejobs::job_single("get_genomewide_SNP_list", memory = "150G",create_shell = TRUE, command = "Rscript get_genomewide_SNP_list.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
