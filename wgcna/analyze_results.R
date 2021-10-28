##
library(sva)
library(lmerTest)
library(SummarizedExperiment)
library(jaffelab)
library(WGCNA)
library(broom)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(readxl)
library(RColorBrewer)
library(sessioninfo)
library(here)

rm(list = ls())


capabilities()

## load data

load("../exprs_cutoff/rse_gene.Rdata", verbose = TRUE)
load("../data/degradation_rse_MDDseq_BiPSeq_BothRegions.Rdata", verbose = TRUE)

## checks
identical(colnames(rse_gene), colnames(cov_rse)) # TRUE

rse_gene <- rse_gene[, rse_gene$PrimaryDx %in% c("Control", "MDD")]
cov_rse <- cov_rse[, cov_rse$PrimaryDx %in% c("Control", "MDD")]

rse_gene$Dx <- droplevels(rse_gene$PrimaryDx)
cov_rse$Dx <- droplevels(rse_gene$PrimaryDx)

rse_gene$Dx<-relevel(rse_gene$Dx, "Control")
rse_gene$Dx <- relevel(rse_gene$Dx, "Control")

table(rse_gene$Dx)

#Control     MDD 
#    387     459 

## add ancestry
#load("../genotype_data/goesHyde_bipolarMdd_Genotypes_n588_mds.rda", verbose = TRUE)

## keep samples with genotypes (remember that rse_gene is summarized experiment, so BrNum is in colData)
#rse_gene <- rse_gene[, rse_gene$BrNum %in% rownames(mds)]
#cov_rse <- cov_rse[, cov_rse$BrNum %in% rownames(mds)]

#### reorders according to rse_gene$BrNum
#mds = mds[rse_gene$BrNum,1:5]

#colData(rse_gene) = cbind(colData(rse_gene), mds)
#colData(cov_rse) = cbind(colData(cov_rse), mds)

###########
#compute RPKM
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')

##### from run_wgcna_coombinedR script this was the model used for wgcna 
colnames(modQsva)

 #[1] "(Intercept)"               "DxMDD"
 #[3] "BrainRegionsACC"           "AgeDeath"
 #[5] "DxMDD:BrainRegionsACC" "SexM"
 #[7] "snpPC1"                    "snpPC2"
 #[9] "snpPC3"                    "mitoRate"
#[11] "rRNA_rate"                 "totalAssignedGene"
#[13] "RIN"                       "ERCCsumLogErr"
#[15] "PC1"                       "PC2"
#[17] "PC3"                       "PC4"
#[19] "PC5"                       "PC6"
#[21] "PC7"                       "PC8"
#[23] "PC9"                       "PC10"
#[25] "PC11"                      "PC12"
#[27] "PC13"                      "PC14"
#[29] "PC15"                      "PC16"
#[31] "PC17"                      "PC18"
#[33] "PC19"                      "PC20"
#[35] "PC21"                      "PC22"
## clean expression
geneExprs = log2(recount::getRPKM(rse_gene, "Length")+1)


load("modQsva.rda", verbose = TRUE)
### regress out after variable 5 (protect 1,2,3,4, 5)
geneExprsClean = cleaningY(geneExprs, modQsva, P=5)

#load(here('wgcna', 'geneExprsClean'), verbose = TRUE)
          
          ##### Load rse data, examine ####

#######end of model used 

##########
## model #
##########

### put back ERCCsumLogErr term
modJoint = model.matrix(~Dx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 +
                          mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr,
                        data=colData(rse_gene))

## counts from degrafation into log2 scale
degExprs = log2(assays(cov_rse)$count+1)
##SVA command
#k = sva::num.sv(degExprs, modJoint)
#unclear why above command gives error
k = num.sv(degExprs, modJoint)

print(k) # 22
### compute principal components from degratation matrix , keep 22 
qSV_mat = prcomp(t(degExprs))$x[,1:k]

colnames(modJoint)
# [1] "(Intercept)"               "DxMDD"
# [3] "BrainRegionsACC"           "AgeDeath"
# [5] "SexM"                      "snpPC1"
# [7] "snpPC2"                    "snpPC3"
# [9] "mitoRate"                  "rRNA_rate"
#[11] "totalAssignedGene"         "RIN"
#[13] "ERCCsumLogErr"             "DxMDD:BrainRegionsACC"

modQsva = cbind(modJoint[,c(1:4,14,5:13)], qSV_mat)

colnames(modQsva)

save(modQsva, file = "modQsva.rda")
load("modQsva.rda", verbose = TRUE)

# [1] "(Intercept)"               "DxMDD"
#  [3] "BrainRegionsACC"           "AgeDeath"
#  [5] "DxMDD:BrainRegionsACC" "SexM"
#  [7] "snpPC1"                    "snpPC2"
#  [9] "snpPC3"                    "mitoRate"
# [11] "rRNA_rate"                 "totalAssignedGene"
# [13] "RIN"                       "ERCCsumLogErr"
# [15] "PC1"                       "PC2"
# [17] "PC3"                       "PC4"
# [19] "PC5"                       "PC6"
# [21] "PC7"                       "PC8"
# [23] "PC9"                       "PC10"
# [25] "PC11"                      "PC12"
# [27] "PC13"                      "PC14"
# [29] "PC15"                      "PC16"
# [31] "PC17"                      "PC18"
# [33] "PC19"                      "PC20"
# [35] "PC21"                      "PC22"


#########################
## load wgcna output ####
#########################

## load
load("rdas/constructed_network_signed_bicor.rda", verbose=TRUE)
# Loading objects:
#   net_list
#   net
#   fNames

#net is the main blockwide module command

### added code to create dendogram from WGCNA Tutorial
mergedColors = labels2colors(net$colors)

pdf("dendrogram_102521.pdf")

plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)

dev.off() 
#########

# get colors LOOK AT WGCNA Instructions
net$colorsLab = labels2colors(net$colors)
colorDat = data.frame(num = net$colors, col = net$colorsLab, 
	stringsAsFactors=FALSE)
colorDat$Label = paste0("ME", colorDat$num)
colorDat = colorDat[order(colorDat$num),]
colorDat = colorDat[!duplicated(colorDat$num),]
colorDat$numGenes = table(net$colors)[as.character(colorDat$num)]

dim(colorDat)
# [1] 17  4
### WILL CHANGE WITH AGE PROTECTED
colorDat

#17 4
#                    num          col Label numGenes
# ENSG00000227232.5    0         grey   ME0    15157
# ENSG00000228794.8    1    turquoise   ME1     2859
# ENSG00000225630.1    2         blue   ME2     1769
# ENSG00000205116.3    3        brown   ME3     1214
# ENSG00000230415.1    4       yellow   ME4      729
# ENSG00000078369.17   5        green   ME5      662
# ENSG00000142609.17   6          red   ME6      600
# ENSG00000188976.10   7        black   ME7      587
# ENSG00000142583.17   8         pink   ME8      579
# ENSG00000196581.10   9      magenta   ME9      338
# ENSG00000162576.16  10       purple  ME10      187
# ENSG00000179546.4   11  greenyellow  ME11      132
# ENSG00000088280.18  12          tan  ME12       94
# ENSG00000060656.19  13       salmon  ME13       94
# ENSG00000187634.11  14         cyan  ME14       85
# ENSG00000077254.14  15 midnightblue  ME15       76
# ENSG00000107404.18  16    lightcyan  ME16       50


######################
## gene ontology
### plit genes by each color 
gList = split(rowData(rse_gene)$EntrezID, 
	factor(net$colorsLab,levels = colorDat$col))
gList = lapply(gList, function(x) as.character(x[!is.na(x)]))
univ = rowData(rse_gene)$EntrezID
univ = as.character(univ[!is.na(univ)])

go = compareCluster(gList, univ = univ,
	OrgDb = "org.Hs.eg.db", ont = "ALL",
	readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 1)
save(go, file = "rdas/go_enrichment_MDD_Control_102511_wgcna.rda")
######
#### split genes with enselble ID 
gList = split(rowData(rse_gene)$ensemblID, 
	factor(net$colorsLab,levels = colorDat$col))
gList = lapply(gList, function(x) as.character(x[!is.na(x)]))
univ = rowData(rse_gene)$ensemblID
univ = as.character(univ[!is.na(univ)])


goDf = as.data.frame(go)
## below specific for MDD classes 
#goCheck = goDf[goDf$Cluster %in% c("red", "pink", "magenta", "royalblue") &
### changed to include nominally significant categories for the association with depression
goCheck = goDf[goDf$Cluster %in% c("lightcyan", "pink") & goDf$qvalue < 0.05,]
goCheck  =goCheck[order(goCheck$pvalue),]
write.csv(goCheck, file = "go_enrichment_MDD_Control_wgcna_CandidateModules.csv")


##############
## associate eigengenes with brain region
m = modQsva[,1:5] # this is what was protected
colnames(m)
#### look up MEs in WGCNA 
MEs = net$MEs
## same order 
colnames(MEs) = colorDat$col[match(colnames(MEs), colorDat$Label)]
MEs = MEs[,colorDat$col] #change order of columns
dim(MEs)
## check## lmer intercept -- random effect -- check what -1 means 
statList = lapply(MEs, function(x) summary(lmer(x ~ m + (1|rse_gene$BrNum) - 1))$coef)

statList[[1]]

#    WILL CHANGE                 Estimate  Std. Error       df    t value
# m(Intercept)               -0.0006301415 0.002246288 641.3682 -0.2805257
# mDxControl                 -0.0080397674 0.003359925 672.6521 -2.3928413
# mBrainRegionsACC            0.0035551921 0.002087894 390.5075  1.7027648
# mDxControl:BrainRegionsACC  0.0093283586 0.003168815 406.2838  2.9438007
#                               Pr(>|t|)
# m(Intercept)               0.779164692
# mDxControl                 0.016991760
# mBrainRegionsACC           0.089407829
# mDxControl:BrainRegionsACC 0.003428175

#updated_092321
#                                 Estimate   Std. Error       df    t value
# m(Intercept)                0.0050057015 4.602188e-03 495.0492  1.0876787
# mDxControl                 -0.0029024699 3.391172e-03 680.2306 -0.8558900
# mBrainRegionsACC           -0.0010713519 2.099795e-03 405.0555 -0.5102174
# mAgeDeath                  -0.0001061467 9.048749e-05 448.5304 -1.1730539
# mDxControl:BrainRegionsACC  0.0081205984 3.137935e-03 415.2745  2.5878795
#                               Pr(>|t|)
# m(Intercept)               0.277266253
# mDxControl                 0.392360007
# mBrainRegionsACC           0.610177080
# mAgeDeath                  0.241396772
# mDxControl:BrainRegionsACC 0.009995445


# modified extract information from lme object
MDDEffect = as.data.frame(t(sapply(statList, function(x) x[2,])))
regionEffect = as.data.frame(t(sapply(statList, function(x) x[3,])))
## now protect age we cannot 
ageEffect = as.data.frame(t(sapply(statList, function(x) x[4,])))  
#interaction effect
intEffect = as.data.frame(t(sapply(statList, function(x) x[5,])))
#colnames(MDDEffect)= colnames(regionEffect) = colnames(ageEffect) = colnames(intEffect) = c(
#	"slope", "se", "df", "t", "pvalue")

colnames(MDDEffect)= colnames(regionEffect) = colnames(ageEffect) = colnames(intEffect) = c(
	"slope", "se", "df", "t", "pvalue")

print_effect <- function(x) { signif(cbind(x, FDR = p.adjust(x$pvalue, 'fdr')), 3) }
print_effect(MDDEffect)


#updated102521
#                  slope      se  df      t   pvalue     FDR
# grey          0.002900 0.00339 680  0.856 0.393000 0.51300
# turquoise     0.001590 0.00167 781  0.948 0.343000 0.48600
# blue          0.002910 0.00174 834  1.670 0.095000 0.21500
# brown         0.000536 0.00332 830  0.161 0.872000 0.91100
# yellow       -0.005440 0.00332 743 -1.640 0.101000 0.21500
# green         0.006430 0.00324 831  1.980 0.047500 0.13500
# red          -0.004450 0.00322 806 -1.380 0.167000 0.28100
# black         0.007840 0.00324 818  2.420 0.015700 0.08880
# pink         -0.011300 0.00333 721 -3.380 0.000761 0.00973
# magenta       0.006820 0.00315 818  2.170 0.030600 0.10900
# purple       -0.001850 0.00330 787 -0.560 0.575000 0.65200
# greenyellow  -0.000341 0.00306 836 -0.111 0.911000 0.91100
# tan           0.006670 0.00310 648  2.150 0.032000 0.10900
# salmon       -0.002490 0.00337 828 -0.739 0.460000 0.55900
# cyan         -0.004310 0.00322 802 -1.340 0.182000 0.28100
# midnightblue -0.004400 0.00285 838 -1.540 0.123000 0.23300
# lightcyan     0.010900 0.00335 774  3.260 0.001140 0.00973



#signif(ageEffect, 3)
print_effect(ageEffect)

#                  slope      se  df      t  pvalue
# grey          0.009330 0.00317 406  2.940 0.00343
# turquoise     0.002460 0.00198 420  1.240 0.21400
# blue          0.003150 0.00234 414  1.340 0.18000
# brown         0.008930 0.00438 424  2.040 0.04210
# yellow        0.001980 0.00392 399  0.505 0.61400
# green         0.004730 0.00391 388  1.210 0.22700
# red           0.001800 0.00349 399  0.516 0.60600
# black         0.001410 0.00427 363  0.329 0.74200
# pink          0.004870 0.00429 404  1.140 0.25700
# magenta       0.006090 0.00398 430  1.530 0.12600
# purple       -0.001000 0.00428 431 -0.235 0.81500
# greenyellow  -0.002170 0.00451 428 -0.480 0.63100
# tan          -0.001710 0.00410 392 -0.417 0.67700
# salmon        0.003180 0.00282 407  1.130 0.26100
# cyan          0.000864 0.00393 435  0.220 0.82600
# midnightblue  0.005490 0.00405 424  1.360 0.17600
# lightcyan     0.004220 0.00430 428  0.979 0.32800
# grey60        0.002750 0.00420 414  0.654 0.51300

#updated092321 Age
#                  slope       se  df      t   pvalue      FDR
# grey         -1.06e-04 9.05e-05 449 -1.170 2.41e-01 2.74e-01
# turquoise     1.12e-04 4.08e-05 435  2.750 6.20e-03 1.17e-02
# blue         -6.52e-05 3.88e-05 429 -1.680 9.35e-02 1.20e-01
# brown         3.56e-04 7.50e-05 432  4.750 2.82e-06 1.60e-05
# yellow        1.38e-04 8.35e-05 417  1.650 9.90e-02 1.20e-01
# green         2.45e-04 7.26e-05 413  3.380 8.01e-04 1.95e-03
# red           1.25e-04 7.56e-05 408  1.660 9.77e-02 1.20e-01
# black        -2.22e-05 7.48e-05 417 -0.297 7.67e-01 8.15e-01
# pink          2.39e-04 8.57e-05 429  2.790 5.46e-03 1.16e-02
# magenta       3.38e-04 7.25e-05 403  4.670 4.19e-06 1.78e-05
# purple        3.60e-04 8.02e-05 444  4.490 9.18e-06 3.12e-05
# greenyellow  -5.04e-04 6.77e-05 427 -7.440 5.50e-13 9.34e-12
# tan          -3.11e-04 8.46e-05 449 -3.670 2.67e-04 7.57e-04
# salmon        1.91e-04 7.64e-05 435  2.500 1.29e-02 2.00e-02
# cyan          5.27e-04 7.67e-05 436  6.870 2.20e-11 1.87e-10
# midnightblue -6.87e-06 6.23e-05 439 -0.110 9.12e-01 9.12e-01
# lightcyan     2.10e-04 8.10e-05 372  2.600 9.71e-03 1.65e-02

print_effect(intEffect)

updated092321
#                  slope      se  df      t pvalue    FDR
# grey          0.008120 0.00314 415  2.590 0.0100 0.0850
# turquoise     0.001240 0.00192 423  0.645 0.5190 0.8020
# blue          0.003820 0.00230 434  1.660 0.0979 0.3420
# brown         0.010100 0.00432 435  2.330 0.0203 0.1150
# yellow        0.011900 0.00359 397  3.310 0.0010 0.0171
# green        -0.001090 0.00424 416 -0.256 0.7980 0.8600
# red           0.002740 0.00395 403  0.694 0.4880 0.8020
# black         0.003720 0.00408 415  0.913 0.3620 0.6830
# pink          0.001380 0.00343 405  0.403 0.6870 0.8600
# magenta      -0.000906 0.00398 402 -0.227 0.8200 0.8600
# purple        0.004170 0.00383 432  1.090 0.2770 0.5880
# greenyellow   0.004790 0.00408 433  1.170 0.2410 0.5860
# tan           0.004140 0.00266 410  1.550 0.1210 0.3420
# salmon       -0.001860 0.00436 437 -0.427 0.6700 0.8600
# cyan          0.001320 0.00388 429  0.341 0.7330 0.8600
# midnightblue -0.000681 0.00385 447 -0.177 0.8600 0.8600
# lightcyan     0.006360 0.00394 361  1.610 0.1080 0.3420
# > 


##################
# make boxplots ##
lab = paste0(substr(rse_gene$BrainRegion,1,4), ":", rse_gene$Dx)
table(lab)
# Amyg:Control     Amyg:MDD sACC:Control     sACC:MDD 
#          187          231          200          228 
### changed BP to MDD labels
lab = factor(lab, levels = c("sACC:Control", "sACC:MDD", "Amyg:Control", "Amyg:MDD"))

pdf("MEs_vs_dx.pdf",w=8,h=6, useDingbats = FALSE)
palette(brewer.pal(4,"Paired"))
par(mar=c(3,6,4,2), cex.axis=2,cex.lab=1.8,cex.main = 1.8)
for(i in 1:ncol(MEs)) {
	boxplot(MEs[,i] ~ lab, outline = FALSE, xlab="",
		ylim = quantile(unlist(MEs),c(0.001,0.999)),main = colnames(MEs)[i],
		names = gsub(":", "\n", levels(lab)), ylab = "Module Eigengene")
	points(MEs[,i] ~ jitter(as.numeric(lab),amount=0.1), pch=21, bg=lab)
	legend("top", c(paste0("Region p=", signif(regionEffect[i,5],3)), 
		paste0("Dx p=", signif(MDDEffect[i,5],3))),cex=1.4)
}
dev.off()

colnames(m)
# [1] "(Intercept)"               "DxControl"                
# [3] "BrainRegionsACC"           "AgeDeath"                 
# [5] "DxControl:BrainRegionsACC"

clean_ME = t(cleaningY(t(MEs), m, P=2))
pdf("clean_MEs_vs_dx.pdf",w=5,h=6, useDingbats = FALSE)
palette(brewer.pal(4,"Paired"))
par(mar=c(3,6,4,2), cex.axis=2,cex.lab=1.8,cex.main = 1.8)
for(i in 1:ncol(clean_ME)) {
	boxplot(clean_ME[,i] ~ rse_gene$Dx, outline = FALSE, xlab="",
		ylim = quantile(unlist(clean_ME),c(0.001,0.999)),main = colnames(MEs)[i],
		ylab = "Module Eigengene (Adj)")
	points(clean_ME[,i] ~ jitter(as.numeric(rse_gene$Dx),amount=0.1), pch=21, bg=lab)
	legend("top", paste0("Dx p=", signif(MDDEffect[i,5],3)),cex=1.4)
}
dev.off()
	

##############################
## enrichment of DEGs #####
#############################
load("../case_control/interaction_model_results.rda")
## remove everything but gene level
rm(outExon_bothRegion, outJxn_bothRegion, outTx_bothRegion)

identical(rownames(outGene_bothRegion), fNames) # TRUE
tt = table(net$colorsLab, outGene_bothRegion$adj.P.Val_dxEffect < 0.05)
tt = tt[rownames(MDDEffect),]

plot(prop.table(tt,1)[,2], -log10(MDDEffect$pvalue))

# manual chi-sq?

deModEnrich = as.data.frame(t(sapply(colorDat$col, function(cc) {
	tab = table(net$colorsLab == cc,outGene_bothRegion$adj.P.Val_dxEffect < 0.05)
	c(chisq.test(tab)$p.value, getOR(tab))
})))
colnames(deModEnrich) = c("Pvalue", "OR")
deModEnrich$numGenes = colorDat$numGenes
deModEnrich$numDE = tt[,2]

plot(-log10(deModEnrich$Pvalue), -log10(MDDEffect$pvalue ))

## do stuff on red

## read in magma
magma = read_excel("../PGC_BP_Magma_table.xlsx", skip = 2)
#updated 
magma = read.table(file = "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/Howard_MDD_MAGMA_genes.csv", sep = ",")


magma = as.data.frame(magma)
magma$FDR_JOINT = p.adjust(magma$P_JOINT, "fdr")
magma$BONF_JOINT = p.adjust(magma$P_JOINT, "bonf")
magma$Module = net$colorsLab[match(magma$GENE, rowData(rse_gene)$EntrezID)]
magma$isSig = factor(ifelse(magma$BONF_JOINT < 0.05, "Yes","No"))
ms = unique(magma$Module)
ms = ms[!is.na(ms)]

## test
magmaEnrich = t(sapply(ms, function(m) {
	tt =table(magma$Module == m, magma$isSig)
	c(getOR(tt), chisq.test(tt)$p.value)
}))
colnames(magmaEnrich) = c("OR", "pvalue")
magmaEnrich = as.data.frame(magmaEnrich)

## add counts
magmaEnrich$MagmaSig = table(magma$Module[magma$BONF_JOINT < 0.05])[rownames(magmaEnrich)]
magmaEnrich$NumGenes = table(magma$Module)[rownames(magmaEnrich)]
write.csv(magmaEnrich, file = "magma_module_enrichment.csv")

magmaSig = magma[magma$BONF_JOINT < 0.05 & magma$GENE %in% rowData(rse_gene)$EntrezID,]

entrezIDs = split(rowData(rse_gene)$EntrezID, net$colorsLab)
entrezOverlap = sapply(entrezIDs, function(x) sum(x %in% magmaSig$GENE,na.rm=TRUE))
entrezPresent = sapply(entrezIDs, function(x) sum(x %in% magma$GENE,na.rm=TRUE))
x = data.frame(numOverlap = entrezOverlap, numPresent = entrezPresent)
x = x[colorDat$col,]
x$numGenes = colorDat$numGenes
x$totGene = nrow(rse_gene)

## pink enrich
mat = matrix(c(5,198, 128, 13894 -5 - 198-128), nr=2)
chisq.test(mat)
getOR(mat)

#magmaSig = magma[which(magma$P_JOINT
## line up
mm = match(rowData(rse_gene)$Symbol, magma$SYMBOL)
magGenes = split(magma$SYMBOL[mm], net$colorsLab)



#####extract genes in module for visualization 
#c("lightcyan", "pink") 


load("rdas/wgcna_signed_TOM-block.1.RData", verbose = TRUE)
TOM.mat = as.matrix(TOM)

#TOM = TOMsimilarityFromExpr(t(geneExprsClean), power = 10)

dim(TOM.mat)    

## clean expression
geneExprs = log2(recount::getRPKM(rse_gene, "Length")+1)

### regress out after variable 5 (protect 1,2,3,4, 5)
geneExprsClean = cleaningY(geneExprs, modQsva, P=5)



# choose module
module = "lightcyan"

module = "pink"

# get list of genes

dim(geneExprsClean)
probes = rownames(geneExprsClean)
#ensembl ids
length(probes)

#25212


inModule = (net$colorsLab==module)
modProbes = probes[inModule]
# Select the corresponding Topological Overlap
modTOM = TOM.mat[inModule, inModule]

dimnames(modTOM) = list(modProbes, modProbes)

kIN=softConnectivity(t(geneExprsClean[modProbes, ]))

selectHubs = (rank (-kIN))

vis = exportNetworkToVisANT(modTOM[selectHubs,selectHubs],
file=paste("VisANTInput-", module, "-top.txt", sep=""),
weighted=TRUE, 
threshold = 0, 
probeToGene= data.frame(rowRanges(rse_gene)$gencodeID, rowRanges(rse_gene)$Symbol) )


#### cytoscape###### 

cyt = exportNetworkToCytoscape(modTOM,
edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
weighted = TRUE,
threshold = 0.02, 
nodeNames = modProbes,
altNodeNames = modGenes)


nodeAttr = moduleColors[inModule])






## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.6.1 Patched (2019-10-31 r77350)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2020-04-06
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version    date       lib source
#  acepack                1.4.1      2016-10-29 [2] CRAN (R 3.6.1)
#  annotate               1.64.0     2019-10-29 [2] Bioconductor
#  AnnotationDbi        * 1.48.0     2019-10-29 [2] Bioconductor
#  askpass                1.1        2019-01-13 [2] CRAN (R 3.6.1)
#  assertthat             0.2.1      2019-03-21 [2] CRAN (R 3.6.1)
#  backports              1.1.5      2019-10-02 [2] CRAN (R 3.6.1)
#  base64enc              0.1-3      2015-07-28 [2] CRAN (R 3.6.1)
#  Biobase              * 2.46.0     2019-10-29 [2] Bioconductor
#  BiocFileCache          1.10.2     2019-11-08 [2] Bioconductor
#  BiocGenerics         * 0.32.0     2019-10-29 [2] Bioconductor
#  BiocManager            1.30.10    2019-11-16 [2] CRAN (R 3.6.1)
#  BiocParallel         * 1.20.1     2019-12-21 [2] Bioconductor
#  biomaRt                2.42.0     2019-10-29 [2] Bioconductor
#  Biostrings             2.54.0     2019-10-29 [2] Bioconductor
#  bit                    1.1-15.2   2020-02-10 [2] CRAN (R 3.6.1)
#  bit64                  0.9-7      2017-05-08 [2] CRAN (R 3.6.1)
#  bitops                 1.0-6      2013-08-17 [2] CRAN (R 3.6.1)
#  blob                   1.2.1      2020-01-20 [2] CRAN (R 3.6.1)
#  boot                   1.3-23     2019-07-05 [3] CRAN (R 3.6.1)
#  broom                * 0.5.4      2020-01-27 [2] CRAN (R 3.6.1)
#  BSgenome               1.54.0     2019-10-29 [2] Bioconductor
#  bumphunter             1.28.0     2019-10-29 [2] Bioconductor
#  cellranger             1.1.0      2016-07-27 [2] CRAN (R 3.6.1)
#  checkmate              2.0.0      2020-02-06 [2] CRAN (R 3.6.1)
#  cli                    2.0.2      2020-02-28 [1] CRAN (R 3.6.1)
#  cluster                2.1.0      2019-06-19 [3] CRAN (R 3.6.1)
#  clusterProfiler      * 3.14.3     2020-01-08 [1] Bioconductor
#  codetools              0.2-16     2018-12-24 [3] CRAN (R 3.6.1)
#  colorout             * 1.2-2      2020-04-02 [1] Github (jalvesaq/colorout@726d681)
#  colorspace             1.4-1      2019-03-18 [2] CRAN (R 3.6.1)
#  cowplot                1.0.0      2019-07-11 [1] CRAN (R 3.6.1)
#  crayon                 1.3.4      2017-09-16 [2] CRAN (R 3.6.1)
#  curl                   4.3        2019-12-02 [2] CRAN (R 3.6.1)
#  data.table             1.12.8     2019-12-09 [2] CRAN (R 3.6.1)
#  DBI                    1.1.0      2019-12-15 [2] CRAN (R 3.6.1)
#  dbplyr                 1.4.2      2019-06-17 [2] CRAN (R 3.6.1)
#  DelayedArray         * 0.12.2     2020-01-06 [2] Bioconductor
#  derfinder              1.20.0     2019-10-29 [2] Bioconductor
#  derfinderHelper        1.20.0     2019-10-29 [2] Bioconductor
#  digest                 0.6.25     2020-02-23 [1] CRAN (R 3.6.1)
#  DO.db                  2.9        2020-04-06 [1] Bioconductor
#  doParallel             1.0.15     2019-08-02 [2] CRAN (R 3.6.1)
#  doRNG                  1.8.2      2020-01-27 [2] CRAN (R 3.6.1)
#  DOSE                   3.12.0     2019-10-29 [1] Bioconductor
#  downloader             0.4        2015-07-09 [2] CRAN (R 3.6.1)
#  dplyr                  0.8.4      2020-01-31 [2] CRAN (R 3.6.1)
#  dynamicTreeCut       * 1.63-1     2016-03-11 [1] CRAN (R 3.6.1)
#  ellipsis               0.3.0      2019-09-20 [2] CRAN (R 3.6.1)
#  enrichplot             1.6.1      2019-12-16 [1] Bioconductor
#  europepmc              0.3        2018-04-20 [1] CRAN (R 3.6.1)
#  fansi                  0.4.1      2020-01-08 [2] CRAN (R 3.6.1)
#  farver                 2.0.3      2020-01-16 [2] CRAN (R 3.6.1)
#  fastcluster          * 1.1.25     2018-06-07 [2] CRAN (R 3.6.1)
#  fastmatch              1.1-0      2017-01-28 [1] CRAN (R 3.6.1)
#  fgsea                  1.12.0     2019-10-29 [1] Bioconductor
#  foreach                1.4.8      2020-02-09 [2] CRAN (R 3.6.1)
#  foreign                0.8-72     2019-08-02 [3] CRAN (R 3.6.1)
#  Formula                1.2-3      2018-05-03 [2] CRAN (R 3.6.1)
#  genefilter           * 1.68.0     2019-10-29 [2] Bioconductor
#  generics               0.0.2      2018-11-29 [2] CRAN (R 3.6.1)
#  GenomeInfoDb         * 1.22.0     2019-10-29 [2] Bioconductor
#  GenomeInfoDbData       1.2.2      2019-10-28 [2] Bioconductor
#  GenomicAlignments      1.22.1     2019-11-12 [2] Bioconductor
#  GenomicFeatures        1.38.2     2020-02-15 [2] Bioconductor
#  GenomicFiles           1.22.0     2019-10-29 [2] Bioconductor
#  GenomicRanges        * 1.38.0     2019-10-29 [2] Bioconductor
#  GEOquery               2.54.1     2019-11-18 [2] Bioconductor
#  ggforce                0.3.1      2019-08-20 [2] CRAN (R 3.6.1)
#  ggplot2                3.2.1      2019-08-10 [2] CRAN (R 3.6.1)
#  ggplotify              0.0.5      2020-03-12 [1] CRAN (R 3.6.1)
#  ggraph                 2.0.1      2020-02-07 [2] CRAN (R 3.6.1)
#  ggrepel                0.8.1      2019-05-07 [2] CRAN (R 3.6.1)
#  ggridges               0.5.2      2020-01-12 [1] CRAN (R 3.6.1)
#  glue                   1.3.2      2020-03-12 [1] CRAN (R 3.6.1)
#  GO.db                  3.10.0     2019-10-28 [2] Bioconductor
#  googledrive            1.0.0      2019-08-19 [1] CRAN (R 3.6.1)
#  GOSemSim               2.12.1     2020-03-19 [1] Bioconductor
#  graphlayouts           0.5.0      2019-08-20 [2] CRAN (R 3.6.1)
#  gridExtra              2.3        2017-09-09 [2] CRAN (R 3.6.1)
#  gridGraphics           0.5-0      2020-02-25 [1] CRAN (R 3.6.1)
#  gtable                 0.3.0      2019-03-25 [2] CRAN (R 3.6.1)
#  Hmisc                  4.3-1      2020-02-07 [2] CRAN (R 3.6.1)
#  hms                    0.5.3      2020-01-08 [2] CRAN (R 3.6.1)
#  htmlTable              1.13.3     2019-12-04 [2] CRAN (R 3.6.1)
#  htmltools              0.4.0      2019-10-04 [2] CRAN (R 3.6.1)
#  htmlwidgets            1.5.1      2019-10-08 [2] CRAN (R 3.6.1)
#  httr                   1.4.1      2019-08-05 [2] CRAN (R 3.6.1)
#  igraph                 1.2.4.2    2019-11-27 [2] CRAN (R 3.6.1)
#  impute                 1.60.0     2019-10-29 [2] Bioconductor
#  IRanges              * 2.20.2     2020-01-13 [2] Bioconductor
#  iterators              1.0.12     2019-07-26 [2] CRAN (R 3.6.1)
#  jaffelab             * 0.99.30    2020-04-02 [1] Github (LieberInstitute/jaffelab@42637ff)
#  jpeg                   0.1-8.1    2019-10-24 [2] CRAN (R 3.6.1)
#  jsonlite               1.6.1      2020-02-02 [2] CRAN (R 3.6.1)
#  knitr                  1.28       2020-02-06 [2] CRAN (R 3.6.1)
#  lattice                0.20-38    2018-11-04 [3] CRAN (R 3.6.1)
#  latticeExtra           0.6-29     2019-12-19 [2] CRAN (R 3.6.1)
#  lazyeval               0.2.2      2019-03-15 [2] CRAN (R 3.6.1)
#  lifecycle              0.2.0      2020-03-06 [1] CRAN (R 3.6.1)
#  limma                  3.42.2     2020-02-03 [2] Bioconductor
#  lme4                 * 1.1-21     2019-03-05 [2] CRAN (R 3.6.1)
#  lmerTest             * 3.1-1      2019-12-13 [1] CRAN (R 3.6.1)
#  locfit                 1.5-9.1    2013-04-20 [2] CRAN (R 3.6.1)
#  magrittr               1.5        2014-11-22 [2] CRAN (R 3.6.1)
#  MASS                   7.3-51.4   2019-03-31 [3] CRAN (R 3.6.1)
#  Matrix               * 1.2-17     2019-03-22 [3] CRAN (R 3.6.1)
#  matrixStats          * 0.55.0     2019-09-07 [2] CRAN (R 3.6.1)
#  memoise                1.1.0      2017-04-21 [2] CRAN (R 3.6.1)
#  mgcv                 * 1.8-30     2019-10-24 [3] CRAN (R 3.6.1)
#  minqa                  1.2.4      2014-10-09 [2] CRAN (R 3.6.1)
#  munsell                0.5.0      2018-06-12 [2] CRAN (R 3.6.1)
#  nlme                 * 3.1-141    2019-08-01 [3] CRAN (R 3.6.1)
#  nloptr                 1.2.1      2018-10-03 [2] CRAN (R 3.6.1)
#  nnet                   7.3-12     2016-02-02 [3] CRAN (R 3.6.1)
#  numDeriv               2016.8-1.1 2019-06-06 [2] CRAN (R 3.6.1)
#  openssl                1.4.1      2019-07-18 [2] CRAN (R 3.6.1)
#  org.Hs.eg.db         * 3.10.0     2019-10-28 [2] Bioconductor
#  pillar                 1.4.3      2019-12-20 [2] CRAN (R 3.6.1)
#  pkgconfig              2.0.3      2019-09-22 [2] CRAN (R 3.6.1)
#  plyr                   1.8.5      2019-12-10 [2] CRAN (R 3.6.1)
#  png                    0.1-7      2013-12-03 [2] CRAN (R 3.6.1)
#  polyclip               1.10-0     2019-03-14 [2] CRAN (R 3.6.1)
#  preprocessCore         1.48.0     2019-10-29 [2] Bioconductor
#  prettyunits            1.1.1      2020-01-24 [2] CRAN (R 3.6.1)
#  progress               1.2.2      2019-05-16 [2] CRAN (R 3.6.1)
#  purrr                  0.3.3      2019-10-18 [2] CRAN (R 3.6.1)
#  qvalue                 2.18.0     2019-10-29 [2] Bioconductor
#  R6                     2.4.1      2019-11-12 [2] CRAN (R 3.6.1)
#  rafalib              * 1.0.0      2015-08-09 [1] CRAN (R 3.6.1)
#  rappdirs               0.3.1      2016-03-28 [2] CRAN (R 3.6.1)
#  RColorBrewer         * 1.1-2      2014-12-07 [2] CRAN (R 3.6.1)
#  Rcpp                   1.0.3      2019-11-08 [2] CRAN (R 3.6.1)
#  RCurl                  1.98-1.1   2020-01-19 [2] CRAN (R 3.6.1)
#  readr                  1.3.1      2018-12-21 [2] CRAN (R 3.6.1)
#  readxl               * 1.3.1      2019-03-13 [2] CRAN (R 3.6.1)
#  recount                1.12.1     2019-11-06 [2] Bioconductor
#  rentrez                1.2.2      2019-05-02 [2] CRAN (R 3.6.1)
#  reshape2               1.4.3      2017-12-11 [2] CRAN (R 3.6.1)
#  rlang                  0.4.5      2020-03-01 [1] CRAN (R 3.6.1)
#  rngtools               1.5        2020-01-23 [2] CRAN (R 3.6.1)
#  rpart                  4.1-15     2019-04-12 [3] CRAN (R 3.6.1)
#  Rsamtools              2.2.2      2020-02-11 [2] Bioconductor
#  RSQLite                2.2.0      2020-01-07 [2] CRAN (R 3.6.1)
#  rstudioapi             0.11       2020-02-07 [2] CRAN (R 3.6.1)
#  rtracklayer            1.46.0     2019-10-29 [2] Bioconductor
#  rvcheck                0.1.8      2020-03-01 [1] CRAN (R 3.6.1)
#  S4Vectors            * 0.24.3     2020-01-18 [2] Bioconductor
#  scales                 1.1.0      2019-11-18 [2] CRAN (R 3.6.1)
#  segmented              1.1-0      2019-12-10 [2] CRAN (R 3.6.1)
#  sessioninfo          * 1.1.1      2018-11-05 [2] CRAN (R 3.6.1)
#  stringi                1.4.6      2020-02-17 [2] CRAN (R 3.6.1)
#  stringr                1.4.0      2019-02-10 [2] CRAN (R 3.6.1)
#  SummarizedExperiment * 1.16.1     2019-12-19 [2] Bioconductor
#  survival               3.1-8      2019-12-03 [2] CRAN (R 3.6.1)
#  sva                  * 3.34.0     2019-10-29 [2] Bioconductor
#  tibble                 3.0.0      2020-03-30 [1] CRAN (R 3.6.1)
#  tidygraph              1.1.2      2019-02-18 [2] CRAN (R 3.6.1)
#  tidyr                  1.0.2      2020-01-24 [2] CRAN (R 3.6.1)
#  tidyselect             1.0.0      2020-01-27 [2] CRAN (R 3.6.1)
#  triebeard              0.3.0      2016-08-04 [1] CRAN (R 3.6.1)
#  tweenr                 1.0.1      2018-12-14 [2] CRAN (R 3.6.1)
#  urltools               1.7.3      2019-04-14 [1] CRAN (R 3.6.1)
#  VariantAnnotation      1.32.0     2019-10-29 [2] Bioconductor
#  vctrs                  0.2.4      2020-03-10 [1] CRAN (R 3.6.1)
#  viridis                0.5.1      2018-03-29 [2] CRAN (R 3.6.1)
#  viridisLite            0.3.0      2018-02-01 [2] CRAN (R 3.6.1)
#  WGCNA                * 1.69       2020-02-28 [1] CRAN (R 3.6.1)
#  withr                  2.1.2      2018-03-15 [2] CRAN (R 3.6.1)
#  xfun                   0.12       2020-01-13 [2] CRAN (R 3.6.1)
#  XML                    3.99-0.3   2020-01-20 [2] CRAN (R 3.6.1)
#  xml2                   1.2.2      2019-08-09 [2] CRAN (R 3.6.1)
#  xtable                 1.8-4      2019-04-21 [2] CRAN (R 3.6.1)
#  XVector                0.26.0     2019-10-29 [2] Bioconductor
#  zlibbioc               1.32.0     2019-10-29 [2] Bioconductor
# 
# [1] /users/fgoes/R/3.6.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/library
# >
