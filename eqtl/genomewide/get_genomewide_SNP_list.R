library(sessioninfo)
library(here)

# create rdas directory
dir.create(here("eqtl", "genomewide", "rdas"), showWarnings = FALSE)
########################

######################## get overlapping snps
## load SNP data

## Amyg / sACC snps
load(here("genotype_data", "goesHyde_bipolarMdd_Genotypes.rda"), verbose = TRUE) # need ~50G ram to load
snpMap$pos_hg19 <- paste0(snpMap$CHR, ":", snpMap$POS)
snpMap_amyg <- snpMap

## DLPFC snps
load("/dcl01/ajaffe/data/lab/brainseq_phase1/genotype_data/brainseq_phase1_Genotypes_n732.rda", verbose = TRUE) # need 60G ram to load
snpMap$pos_hg19 <- paste0(snpMap$CHR, ":", snpMap$POS)
rm(snp, mds)

table(snpMap_amyg$rsNumGuess %in% snpMap$rsNumGuess)
# FALSE    TRUE 
# 342276 8549784 
table(!is.na(snpMap_amyg$chr_hg38))
snpMap_amyg <- snpMap_amyg[!is.na(snpMap_amyg$chr_hg38),]
snpMap <- snpMap[!is.na(snpMap$chr_hg38),]

snpMap_amyg$SNP_hg38 <- paste0(snpMap_amyg$chr_hg38, ":", snpMap_amyg$pos_hg38, ":", snpMap_amyg$COUNTED, ":", snpMap_amyg$ALT)
snpMap$SNP_hg38 <- paste0(snpMap$chr_hg38, ":", snpMap$pos_hg38, ":", snpMap$COUNTED, ":", snpMap$ALT)

common_snp_hg38 <- intersect(snpMap_amyg$SNP_hg38, snpMap$SNP_hg38)
length(common_snp_hg38)
# [1] 49852

#8M
nrow(snpMap_amyg) - length(unique(snpMap_amyg$SNP_hg38)) #[1] 3
#6M
nrow(snpMap) - length(unique(snpMap$SNP_hg38)) #[1] 59

snpMap_amyg_match <- snpMap_amyg[snpMap_amyg$SNP_hg38 %in% snpMap$SNP_hg38,]
nrow(snpMap_amyg_match)
snpMap_match <- snpMap[match(snpMap_amyg_match$SNP_hg38,snpMap$SNP_hg38),]
nrow(snpMap_match)
## pretty sure these are the same!(except NAs)
table(snpMap_amyg_match$rsNumGuess == snpMap_match$rsNumGuess)

## whats up w/ rsNumGuess?
table(is.na(snpMap_amyg_match$rsNumGuess))
table(is.na(snpMap_match$rsNumGuess))

common_rsNumGuess <- intersect(unique(snpMap_amyg$rsNumGuess), unique(snpMap$name))
length(common_rsNumGuess)
# [1] 67628

snpMap_mismatch <- snpMap_amyg[
  (snpMap_amyg$rsNumGuess %in% common_rsNumGuess) & 
    !(snpMap_amyg$SNP_hg38 %in% common_snp_hg38) & 
    !is.na(snpMap_amyg$rsNumGuess),]

nrow(snpMap_mismatch)
head(snpMap_mismatch)

message(Sys.time(), " Get snp index")
# ind <- which(snpMap$SNP %in% snpMap_amyg$SNP)
ind <- which(snpMap$SNP_hg38 %in% common_snp_hg38)
message(Sys.time(), " Filter snpMap")
snpMap <- snpMap[ind, ]

message(Sys.time(), " Save snpMapKeep data")
snpMapKeep <- snpMap
save(snpMapKeep, file = here("eqtl", "genomewide", "rdas", "overlappingSNPs.rda"))


message(Sys.time(), " Get amyg snp index")
indAmyg <-  match(snpMap$SNP, snpMap_amyg$SNP)
message(Sys.time(), " Filter snpMap_amyg")
snpMap_amyg <- snpMap_amyg[indAmyg, ]


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ## check
message(Sys.time(), " Run identical check POS")
stopifnot(identical(snpMap$POS, snpMap_amyg$POS))
message(Sys.time(), " Run identical check newRef")
stopifnot(identical(snpMap$newRef, snpMap_amyg$newRef))
message(Sys.time(), " Run identical check rsNumGuess")
stopifnot(identical(snpMap$rsNumGuess, snpMap_amyg$rsNumGuess))


# sgejobs::job_single("get_genomewide_SNP_list", memory = "90G",create_shell = TRUE, command = "Rscript get_genomewide_SNP_list.R")


