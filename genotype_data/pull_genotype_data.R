########################

## read packages
library(jaffelab)
library(readr)
library(SummarizedExperiment)
library(stringr)
library(GenomicRanges)
library(sessioninfo)
library(SNPlocs.Hsapiens.dbSNP149.GRCh38)
library(here)

## load data
load(here("data","rse_gene_GoesZandi.rda"), verbose = TRUE)
pd = colData(rse_gene)
pd_genoSamples <- unique(pd[,c("BrNum", "genoSample")])
rownames(pd_genoSamples) <- NULL

################
## read in #####
message("\n ***** Read in genotypes *****")
## read in genotypes
newbfile = "goesHyde_mdd_Genotypes_maf01_geno10_hwe1e6"
genotypes  = read_delim(here("genotype_data",paste0(newbfile, ".traw")), delim="\t")


snp = as.data.frame(genotypes[,-(1:6)])
snp_BrNum = pd_genoSamples$BrNum[match(colnames(snp), pd_genoSamples$genoSample)]
length(unique(snp_BrNum)) == nrow(pd_genoSamples)

colnames(snp) = snp_BrNum

# make map
snpMap = as.data.frame(genotypes[,(1:6)])
snpMap$CHR[snpMap$CHR=="23"] = "X"
rm(genotypes)
###########################
# fix SNPs that are in/dels
message("\n***** Fix genotypes ***** ", Sys.time())
ncRef= nchar(snpMap$ALT)
ncCount= nchar(snpMap$COUNTED)
snpMap$Type = "SNV"
snpMap$Type[ncRef > ncCount] = "Deletion"
snpMap$Type[ncRef < ncCount] = "Insertion"

snpMap$newRef = snpMap$ALT
snpMap$newCount = snpMap$COUNTED

## update variant types
# deletion
dIndex = which(snpMap$Type=="Deletion")
snpMap$newRef[dIndex] = sapply(dIndex, function(i) 
  sub(snpMap$newCount[i], "",  snpMap$newRef[i]))
snpMap$newCount[dIndex] = "-"

# insertion
iIndex = which(snpMap$Type=="Insertion")
snpMap$newCount[iIndex] = sapply(iIndex, function(i)
  sub(snpMap$newRef[i], "", snpMap$newCount[i]))
snpMap$newRef[iIndex] = "-"
head(snpMap[snpMap$Type != "SNV",],10)

rm(dIndex, iIndex,ncRef, ncCount)
#### 
## dbsnp 142 for rs number
message("\n***** dbsnp 142 ***** ",Sys.time())
rs = read_delim("/dcs01/ajaffe/Annotation/dbsnp142_common.txt", delim = "\t")
rs$class = factor(rs$class, c("single", "deletion", "insertion",
                              "in-del", "microsatellite", "mnp"))
## put in SNV order
message("put in SNV order")
rs = rs[order(rs$class, rs$chromStart),]

## try matching
message("try matching")
snpMap$name = NA
snpMap$name[snpMap$Type != "SNV"] = rs$name[match(
  paste0("chr", snpMap$CHR, ":", snpMap$POS)[snpMap$Type != "SNV"], 
  paste0(rs$chrom, ":", rs$chromStart))]
snpMap$name[snpMap$Type == "SNV"] = rs$name[match(
  paste0("chr", snpMap$CHR, ":", snpMap$POS)[snpMap$Type == "SNV"], 
  paste0(rs$chrom, ":", rs$chromEnd))]

rm(rs)

load("dbSNP.Rdata", verbose = TRUE)
dbSnp149 = dbSnp149[dbSnp149$RefSNP_id %in% snpMap$name]

### match to SNPs
message("match to SNPs")
snpMapGR = GRanges(paste0("chr", snpMap$CHR), 
                   IRanges(snpMap$POS,width=1,names=snpMap$SNP))
oo = findOverlaps(snpMapGR, dbSnp142)

snpMap$rsNumGuess = NA
snpMap$rsNumGuess[queryHits(oo)] = dbSnp142$RefSNP_id[subjectHits(oo)]
rm(dbSnp142)
## if name doesnt exist
message("if name doesnt exist")
snpMap$name[is.na(snpMap$name)] = snpMap$rsNumGuess[is.na(snpMap$name)]

## lastly add existing rs numbers
message("add existing rs numbers")
snpMap$name[grepl("^rs", snpMap$SNP)] = ss(snpMap$SNP[grepl("^rs", snpMap$SNP)], ":")

# make SNP id the name for those still missing
snpMap$name[is.na(snpMap$name)] = snpMap$SNP[is.na(snpMap$name)]


## add coordinates
mm = match(snpMap$name, dbSnp149$RefSNP_id)
snpMap$chr_hg38 = as.character(seqnames(dbSnp149))[mm]
snpMap$pos_hg38 = start(dbSnp149)[mm]

## lift
bimBed = GRanges(paste0("chr", snpMap$CHR), 
                 IRanges(snpMap$POS, width=1,names=snpMap$SNP))
bimBed = bimBed[names(bimBed) %in% snpMap$SNP[is.na(snpMap$chr_hg38)]]
chain = import.chain("/dcl01/lieber/ajaffe/Brain/Imputation/Merged/hg19ToHg38.over.chain")
lifted = unlist(GRangesList(liftOver(bimBed,chain)))

## rematch
mm2 = match(snpMap$SNP, names(lifted))
snpMap$chr_hg38[!is.na(mm2)] = as.character(seqnames(lifted))[mm2[!is.na(mm2)]]
snpMap$pos_hg38[!is.na(mm2)] = start(lifted)[mm2[!is.na(mm2)]]

rm(mm, mm2)
#### read in MDS
mds = read.table(paste0(newbfile, ".mds"), 
                 header=TRUE,as.is=TRUE)
mds$genoSample = paste0(mds$FID,"_",mds$IID)
message("All genoSamples in mds: ", all(pd_genoSamples$genoSample %in%  mds$genoSample))
mds$BrNum <- pd_genoSamples$BrNum[match(mds$genoSample, pd_genoSamples$genoSample)]
rownames(mds) = mds$BrNum

mds = mds[,(4:13)]
colnames(mds) = paste0("snpPC",1:ncol(mds))


#############
## save #####
message("\n***** Save data ***** ", Sys.time())
save(mds, snp, snpMap, compress=TRUE,
     file = "goesHyde_bipolarMdd_Genotypes.rda")
save(mds, compress=TRUE,
     file = "goesHyde_bipolarMdd_Genotypes_mds.rda")

# sgejobs::job_single('pull_genotype_data', create_shell = TRUE, queue= 'bluejay', memory = '150G', command = "Rscript pull_genotype_data.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
