########################

## read packages
library(jaffelab)
library(readr)
library(SummarizedExperiment)
library(stringr)
library(GenomicRanges)

## load data
load("../exprs_cutoff/rse_gene.Rdata")
pd = colData(rse_gene)

## get BrNums
BrNums = BrNumOriginal = pd$BrNum
BrUniqueOriginal = unique(pd$BrNum)

### temporarily mislabel BrNumbers to pull correct genotype info
# - swap DNA BrNum labels from Br2260 and Br2473 when pulling genotype data
BrNums[which(BrNums %in% c("Br2260","Br2473"))]
BrNums[which(BrNums %in% c("Br2260","Br2473"))] = c("Br2473","Br2260","Br2473","Br2260")
BrNums[which(BrNums %in% c("Br2260","Br2473"))]
# - swap DNA BrNum labels from Br2301 and Br2538 when pulling genotype data
BrNums[which(BrNums %in% c("Br2301","Br2538"))]
BrNums[which(BrNums %in% c("Br2301","Br2538"))] = c("Br2538","Br2301","Br2301")
BrNums[which(BrNums %in% c("Br2301","Br2538"))]
# - change DNA label from Br2385 to Br2533 when pulling genotypes
BrNums[which(BrNums %in% c("Br2385","Br2533"))]
BrNums[which(BrNums %in% c("Br2385","Br2533"))] = c("Br2385","Br2385")
BrNums[which(BrNums %in% c("Br2385","Br2533"))]
# - swap DNA BrNum labels from Br5434 and Br5435 when pulling genotype data
BrNums[which(BrNums %in% c("Br5434","Br5435"))]
BrNums[which(BrNums %in% c("Br5434","Br5435"))] = c("Br5435","Br5434")
BrNums[which(BrNums %in% c("Br5434","Br5435"))]

pd$BrNum = BrNums
BrUniqueSwapped = unique(BrNums)

############# 
# fam file ##

### read in fam
bfile="/dcl01/lieber/ajaffe/Brain/Imputation/Merged/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_QuadsPlus"
fam = read.table(paste0(bfile, ".fam"), as.is=TRUE)
colnames(fam) = c("FID", "IID", "MID", "PID", "SEX","PHENO")
## fix 5 digit Brain numbers
## fix 5 digit Brain numbers
fam$BrNum = ss(fam$FID,"_")
fam$Batch = ss(fam$FID,"_",2)
fam$BrNum[nchar(fam$BrNum)==7] = paste0("Br", 
	substr(fam$BrNum[nchar(fam$BrNum)==7], 4, 7))


table(unique(pd$BrNum) %in% fam$BrNum) # keep samples w/ genotypes
unique(pd$BrNum) [ ! unique(pd$BrNum) %in% fam$BrNum] # keep samples w/ genotypes
fam = fam[!duplicated(fam$BrNum),] # remove 650 if they have 1M

famOut = fam[which(fam$BrNum %in% pd$BrNum),]
write.table(famOut[,1:2], "samples_to_extract.txt",
	col.names=FALSE, row.names=FALSE, quote=FALSE)
	
#### overall extraction
# bfile = "/dcs01/ajaffe/Imputation/Merged/LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_maf005_hwe1e6_geno10"
newbfile = "goesHyde_mdd_Genotypes_n593_maf01_geno10_hwe1e6"

## extract
system(paste("plink --bfile", bfile, 
	"--keep samples_to_extract.txt --geno 0.1 --maf 0.01 --hwe 0.000001 --make-bed --out", newbfile))

# ## independent and cluster
system(paste("plink --bfile", newbfile, "--maf 0.1 --indep 100 10 1.25 --out", newbfile))

## MDS components	
system(paste0("plink --bfile ", newbfile, 
	" --cluster --mds-plot 10 --extract ",newbfile, ".prune.in --out ", newbfile))

# ## A transpose
system(paste("plink --bfile", newbfile,
	"--recode A-transpose --out", newbfile))
	
################
## read in #####

## read in genotypes
genotypes  = read_delim(paste0(newbfile, ".traw"), delim="\t")

snp = as.data.frame(genotypes[,-(1:6)])
colnames(snp) = ifelse(grepl("^Br", ss(colnames(snp), "_")),
			ss(colnames(snp), "_"), ss(colnames(snp), "_",2))

# make map
snpMap = as.data.frame(genotypes[,(1:6)])
snpMap$CHR[snpMap$CHR=="23"] = "X"

###########################
# fix SNPs that are in/dels 
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

#### 
## dbsnp 142 for rs number
rs = read_delim("/dcs01/ajaffe/Annotation/dbsnp142_common.txt", delim = "\t")
rs$class = factor(rs$class, c("single", "deletion", "insertion",
	"in-del", "microsatellite", "mnp"))
## put in SNV order
rs = rs[order(rs$class, rs$chromStart),]

## try matching
snpMap$name = NA
snpMap$name[snpMap$Type != "SNV"] = rs$name[match(
	paste0("chr", snpMap$CHR, ":", snpMap$POS)[snpMap$Type != "SNV"], 
	paste0(rs$chrom, ":", rs$chromStart))]
snpMap$name[snpMap$Type == "SNV"] = rs$name[match(
	paste0("chr", snpMap$CHR, ":", snpMap$POS)[snpMap$Type == "SNV"], 
	paste0(rs$chrom, ":", rs$chromEnd))]

### lets start w/ hg19 for the remaining ~3M
library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
dbSnp142_list = mclapply(paste0("ch",c(1:22,"X")), function(x) {
	cat(".")
	y = getSNPlocs(x, as.GRanges=TRUE)
	seqlevels(y) = gsub("ch", "chr", seqlevels(y))
	y$RefSNP_id = paste0("rs", y$RefSNP_id)
	return(y)
},mc.cores=6)
dbSnp142 = unlist(GRangesList(dbSnp142_list))

### match to SNPs
snpMapGR = GRanges(paste0("chr", snpMap$CHR), 
	IRanges(snpMap$POS,width=1,names=snpMap$SNP))
oo = findOverlaps(snpMapGR, dbSnp142)

snpMap$rsNumGuess = NA
snpMap$rsNumGuess[queryHits(oo)] = dbSnp142$RefSNP_id[subjectHits(oo)]

## if name doesnt exist
snpMap$name[is.na(snpMap$name)] = snpMap$rsNumGuess[is.na(snpMap$name)]

## lastly add existing rs numbers
snpMap$name[grepl("^rs", snpMap$SNP)] = ss(snpMap$SNP[grepl("^rs", snpMap$SNP)], ":")

# make SNP id the name for those still missing
snpMap$name[is.na(snpMap$name)] = snpMap$SNP[is.na(snpMap$name)]

#####################################################
############# hg38 coordinates ######################

library(SNPlocs.Hsapiens.dbSNP149.GRCh38)
dbSnp149_list = mclapply(c(1:22,"X"), function(x) {
	cat(".")
	y = snpsBySeqname(SNPlocs.Hsapiens.dbSNP149.GRCh38, x)
	seqlevels(y) = paste0("chr", seqlevels(y))
	y = y[y$RefSNP_id %in% snpMap$name]
	return(y)
},mc.cores=6)

dbSnp149 = unlist(GRangesList(lapply(dbSnp149_list, as, "GRanges")))

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


#### read in MDS
mds = read.table(paste0(newbfile, ".mds"), 
	header=TRUE,as.is=TRUE)
mds$BrNum = ss(mds$FID,"_")
mds$BrNum[nchar(mds$BrNum)==7] = paste0("Br", 
	substr(mds$BrNum[nchar(mds$BrNum)==7], 4, 7))
rownames(mds) = mds$BrNum

mds = mds[,(4:13)]
colnames(mds) = paste0("snpPC",1:ncol(mds))


##########################
## correct BrNumbers #####

## confirm order stayed the same throughout
identical(rownames(mds), colnames(snp))
identical(rownames(mds), BrUniqueSwapped)

## reset to original correct BrNums
mm_brain = match( colnames(snp),BrUniqueSwapped)
rownames(mds) = colnames(snp) = BrUniqueOriginal[mm_brain]

#############
## save #####
save(mds, snp, snpMap, compress=TRUE,
	file = "goesHyde_bipolarMdd_Genotypes_n593.rda")
save(mds, compress=TRUE,
	file = "goesHyde_bipolarMdd_Genotypes_n593_mds.rda")
	
	
