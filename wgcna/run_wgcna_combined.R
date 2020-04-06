#######
library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(WGCNA)
library(sva)
library(sessioninfo)


## multithread
allowWGCNAThreads(8)

dir.create("rdas", showWarnings = FALSE)

## load data
load("../exprs_cutoff/rse_gene.Rdata", verbose = TRUE)
load("../data/degradation_rse_MDDseq_BiPSeq_BothRegions.Rdata", verbose = TRUE)

## checks
identical(colnames(rse_gene), colnames(cov_rse)) # TRUE
table(rse_gene$PrimaryDx)
# MDD Control Bipolar
# 463     385     245

table(cov_rse$PrimaryDx)

rse_gene <- rse_gene[, rse_gene$PrimaryDx %in% c("Control", "MDD")]
cov_rse <- cov_rse[, cov_rse$PrimaryDx %in% c("Control", "MDD")]
table(rse_gene$PrimaryDx)
# MDD Control Bipolar
# 463     385       0

rse_gene$Dx <- droplevels(rse_gene$PrimaryDx)
cov_rse$Dx <- droplevels(cov_rse$PrimaryDx)

table(rse_gene$Dx)
# MDD Control
# 463     385

## add ancestry
load("../genotype_data/goesHyde_bipolarMdd_Genotypes_n588_mds.rda", verbose = TRUE)
class(mds)
dim(mds)
corner(mds)
table(rse_gene$BrNum %in% rownames(mds))

# FALSE  TRUE
# 12   836
## 6 individual, 12 brain regions absent in mds file

table(rse_gene$BrNum %in% rownames(mds), rse_gene$BrainRegion)
# Amygdala sACC
# FALSE        6    6
# TRUE       415  421

## keep samples with genotypes
rse_gene <- rse_gene[, rse_gene$BrNum %in% rownames(mds)]
cov_rse <- cov_rse[, cov_rse$BrNum %in% rownames(mds)]

dim(rse_gene)
# [1] 25212   836
addmargins(table(rse_gene$Dx, rse_gene$BrainRegion))
# Amygdala sACC Sum
# MDD          234  227 461
# Control      181  194 375
# Sum          415  421 836

#### reorders according to rse_gene$BrNum
mds = mds[rse_gene$BrNum,1:5]
dim(mds)
# [1] 836   5

#colnames(mds) = paste0("snpPC", 1:5)
colData(rse_gene) = cbind(colData(rse_gene), mds)
colData(cov_rse) = cbind(colData(cov_rse), mds)

###########
#compute RPKM
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')


##########
## model #
##########

### came from qSV_model_DE_analysis.R
# modJoint = model.matrix(~PrimaryDx + AgeDeath + Sex + mitoRate + rRNA_rate +
                            # totalAssignedGene + RIN, data = colData(rse_gene))


### removed ERCCsumLogErr term
modJoint = model.matrix(~Dx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 +
	mitoRate + rRNA_rate + totalAssignedGene + RIN,
	data=colData(rse_gene))


### counts from degrafation into log2 scale
degExprs = log2(assays(cov_rse)$count+1)
k = num.sv(degExprs, modJoint)
print(k) # 22
qSV_mat = prcomp(t(degExprs))$x[,1:k]

## join and move around region, dx and interaction for cleaning

colnames(modJoint)
# [1] "(Intercept)"               "DxControl"
# [3] "BrainRegionsACC"           "AgeDeath"
# [5] "SexM"                      "snpPC1"
# [7] "snpPC2"                    "snpPC3"
# [9] "mitoRate"                  "rRNA_rate"
# [11] "totalAssignedGene"         "RIN"
# [13] "DxControl:BrainRegionsACC"

modQsva = cbind(modJoint[,c(1:3,13,4:12)], qSV_mat)

colnames(modQsva)
# [1] "(Intercept)"               "DxControl"
# [3] "BrainRegionsACC"           "DxControl:BrainRegionsACC"
# [5] "AgeDeath"                  "SexM"
# [7] "snpPC1"                    "snpPC2"
# [9] "snpPC3"                    "mitoRate"
# [11] "rRNA_rate"                 "totalAssignedGene"
# [13] "RIN"                       "PC1"
# [15] "PC2"                       "PC3"
# [17] "PC4"                       "PC5"
# [19] "PC6"                       "PC7"
# [21] "PC8"                       "PC9"
# [23] "PC10"                      "PC11"
# [25] "PC12"                      "PC13"
# [27] "PC14"                      "PC15"
# [29] "PC16"                      "PC17"
# [31] "PC18"                      "PC19"
# [33] "PC20"                      "PC21"
# [35] "PC22"


## clean expression
geneExprs = log2(recount::getRPKM(rse_gene, "Length")+1)

### regress out after variable 4 (protect 1,2,3,4)
geneExprsClean = cleaningY(geneExprs, modQsva, P=4)

#########################
## get power
powers <- c(1:10, seq(from = 12, to=20, by=2))
sftthresh1 <- pickSoftThreshold(t(geneExprsClean), powerVector = powers,
                               networkType = "signed", verbose = 5)
cat(sftthresh1$powerEstimate)
save(sftthresh1, file = "rdas/power_object.rda")

## run wgcna
net = blockwiseModules(t(geneExprsClean), power = sftthresh1$powerEstimate,
                            networkType = "signed", minModuleSize = 30,corType="bicor",
                            reassignThreshold = 0, mergeCutHeight = 0.25,
                            numericLabels = TRUE, pamRespectsDendro = FALSE,
                            saveTOMs = TRUE, verbose = 5, maxBlockSize = 30000,
                            saveTOMFileBase = "rdas/wgcna_signed_TOM")
fNames = rownames(geneExprs)
save(net, fNames, file = "rdas/constructed_network_signed_bicor.rda")

########################
## by region remove ERCCsumLogErr
modRegion =  model.matrix(~Dx + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 +
	mitoRate + rRNA_rate + totalAssignedGene + RIN,
	data=colData(rse_gene))

colnames(modRegion)
# [1] "(Intercept)"       "DxControl"         "AgeDeath"
# [4] "SexM"              "snpPC1"            "snpPC2"
# [7] "snpPC3"            "mitoRate"          "rRNA_rate"
# [10] "totalAssignedGene" "RIN"

###jaffe lab function
rIndexes = splitit(rse_gene$BrainRegion)

## clean by region ## changed P=3 to P=2 since not clear why to protect AgeDeath
geneExprs_list = mclapply(rIndexes, function(ii) {
	degExprs = log2(assays(cov_rse[,ii])$count+1)
	k = num.sv(degExprs, modRegion[ii,])
	m = cbind(modRegion[ii,], prcomp(t(degExprs))$x[,1:k])
	cleaningY(geneExprs[,ii], m, P=2)
},mc.cores=2)

## threshold
thresh_list = mclapply(geneExprs_list, function(y) {
	pickSoftThreshold(t(y), powerVector = powers,
                               networkType = "signed", verbose = 5)
},mc.cores=2)

## networks
net_list = lapply(1:2, function(i) {
	blockwiseModules(t(geneExprs_list[[i]]), power = thresh_list[[i]]$powerEstimate,
                            networkType = "signed", minModuleSize = 30,corType="bicor",
                            reassignThreshold = 0, mergeCutHeight = 0.25,
                            numericLabels = TRUE, pamRespectsDendro = FALSE,
                            saveTOMs = TRUE, verbose = 5, maxBlockSize = 30000,
                            saveTOMFileBase = paste0("rdas/wgcna_signed_TOM_region",
								names(rIndexes)[i]))
})
save(net_list, net, fNames, file = "rdas/constructed_network_signed_bicor.rda")


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

