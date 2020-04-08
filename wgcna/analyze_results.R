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


## load data

load("../exprs_cutoff/rse_gene.Rdata", verbose = TRUE)
load("../data/degradation_rse_MDDseq_BiPSeq_BothRegions.Rdata", verbose = TRUE)

## checks
identical(colnames(rse_gene), colnames(cov_rse)) # TRUE

rse_gene <- rse_gene[, rse_gene$PrimaryDx %in% c("Control", "MDD")]
cov_rse <- cov_rse[, cov_rse$PrimaryDx %in% c("Control", "MDD")]

rse_gene$Dx <- droplevels(rse_gene$PrimaryDx)
cov_rse$Dx <- droplevels(cov_rse$PrimaryDx)

## add ancestry
load("../genotype_data/goesHyde_bipolarMdd_Genotypes_n588_mds.rda", verbose = TRUE)

## keep samples with genotypes (remember that rse_gene is summarized experiment, so BrNum is in colData)
rse_gene <- rse_gene[, rse_gene$BrNum %in% rownames(mds)]
cov_rse <- cov_rse[, cov_rse$BrNum %in% rownames(mds)]

#### reorders according to rse_gene$BrNum
mds = mds[rse_gene$BrNum,1:5]

colData(rse_gene) = cbind(colData(rse_gene), mds)
colData(cov_rse) = cbind(colData(cov_rse), mds)

###########
#compute RPKM
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')

##########
## model #
##########

### removed ERCCsumLogErr term
modJoint = model.matrix(~Dx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 +
                          mitoRate + rRNA_rate + totalAssignedGene + RIN,
                        data=colData(rse_gene))

## counts from degrafation into log2 scale
degExprs = log2(assays(cov_rse)$count+1)
##SVA command
k = sva::num.sv(degExprs, modJoint)
print(k) # 22
### compute principal components from degratation matrix , keep 22 
qSV_mat = prcomp(t(degExprs))$x[,1:k]

colnames(modJoint)
# [1] "(Intercept)"               "DxControl"
# [3] "BrainRegionsACC"           "AgeDeath"
# [5] "SexM"                      "snpPC1"
# [7] "snpPC2"                    "snpPC3"
# [9] "mitoRate"                  "rRNA_rate"
# [11] "totalAssignedGene"         "RIN"
# [13] "DxControl:BrainRegionsACC"

modQsva = cbind(modJoint[,c(1:4,13,5:12)], qSV_mat)

colnames(modQsva)

#  [1] "(Intercept)"               "DxControl"
#  [3] "BrainRegionsACC"           "AgeDeath"
#  [5] "DxControl:BrainRegionsACC" "SexM"
#  [7] "snpPC1"                    "snpPC2"
#  [9] "snpPC3"                    "mitoRate"
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
# >
#########################
## load wgcna output ####
#########################

## load
load("rdas/constructed_network_signed_bicor.rda", verbose=TRUE)



### added code to create dendogram from WGCNA Tutorial
mergedColors = labels2colors(net$colors)

plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)

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
# [1] 18  4
### WILL CHANGE WITH AGE PROTECTED
colorDat

#                    num          col Label numGenes
# ENSG00000227232.5    0         grey   ME0    15381
# ENSG00000228794.8    1    turquoise   ME1     2777
# ENSG00000225630.1    2         blue   ME2     1720
# ENSG00000162572.19   3        brown   ME3     1626
# ENSG00000176022.4    4       yellow   ME4      820
# ENSG00000142609.17   5        green   ME5      686
# ENSG00000142583.17   6          red   ME6      564
# ENSG00000107404.18   7        black   ME7      310
# ENSG00000196581.10   8         pink   ME8      222
# ENSG00000162576.16   9      magenta   ME9      199
# ENSG00000070886.11  10       purple  ME10      192
# ENSG00000162545.5   11  greenyellow  ME11      159
# ENSG00000158286.12  12          tan  ME12      134
# ENSG00000088280.18  13       salmon  ME13      108
# ENSG00000077254.14  14         cyan  ME14       97
# ENSG00000187634.11  15 midnightblue  ME15       97
# ENSG00000179546.4   16    lightcyan  ME16       64
# ENSG00000117155.16  17       grey60  ME17       56
# >

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
save(go, file = "rdas/go_enrichment_MDD_Control_wgcna.rda")


goDf = as.data.frame(go)
## below specific for bipolar classes 
#goCheck = goDf[goDf$Cluster %in% c("red", "pink", "magenta", "royalblue") &
### changed to include nominally significant categories for the association with depression
goCheck = goDf[goDf$Cluster %in% c("grey60", "yellow", "black", "pink", "salmon", "blue", "grey", "cyan", "turquoise", "red") &
	goDf$qvalue < 0.05,]
goCheck  =goCheck[order(goCheck$pvalue),]
write.csv(goCheck, file = "go_enrichment_MDD_Control_wgcna_three_TEST_CandidateModules.csv")
##############
## associate eigengenes with brain region
m = modQsva[,1:5] # this is what was protected

#### look up MEs in WGCNA 
MEs = net$MEs
## same order 
colnames(MEs) = colorDat$col[match(colnames(MEs), colorDat$Label)]
MEs = MEs[,colorDat$col]

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

# modified extract information from lme object
MDDEffect = as.data.frame(t(sapply(statList, function(x) x[2,])))
regionEffect = as.data.frame(t(sapply(statList, function(x) x[3,])))
## now protect age we cannot 
ageEffect = as.data.frame(t(sapply(statList, function(x) x[4,])))  
#interaction effect
intEffect = as.data.frame(t(sapply(statList, function(x) x[5,])))
#colnames(MDDEffect)= colnames(regionEffect) = colnames(ageEffect) = colnames(intEffect) = c(
#	"slope", "se", "df", "t", "pvalue")

colnames(MDDEffect)= colnames(regionEffect) = colnames(intEffect) = c(
	"slope", "se", "df", "t", "pvalue")

print_effect <- function(x) { signif(cbind(x, FDR = p.adjust(x$pvalue, 'fdr')), 3) }
print_effect(MDDEffect)

#                 slope      se  df       t   pvalue
# grey         -0.00804 0.00336 673 -2.3900 0.017000
# turquoise    -0.00360 0.00164 790 -2.1900 0.028500
# blue         -0.00429 0.00174 827 -2.4600 0.013900
# brown         0.00206 0.00338 819  0.6100 0.542000
# yellow       -0.01030 0.00309 809 -3.3500 0.000854
# green        -0.00521 0.00312 802 -1.6700 0.095400
# red           0.00700 0.00336 712  2.0900 0.037400
# black        -0.01010 0.00341 799 -2.9700 0.003050
# pink         -0.00956 0.00333 815 -2.8700 0.004170
# magenta      -0.00027 0.00329 793 -0.0821 0.935000
# purple        0.00563 0.00329 820  1.7100 0.087700
# greenyellow  -0.00516 0.00327 831 -1.5800 0.116000
# tan          -0.00504 0.00338 788 -1.4900 0.137000
# salmon       -0.00834 0.00312 657 -2.6700 0.007740
# cyan          0.00681 0.00290 829  2.3500 0.019200
# midnightblue  0.00150 0.00333 795  0.4520 0.652000
# lightcyan    -0.00281 0.00316 830 -0.8910 0.373000
# grey60       -0.01090 0.00319 823 -3.4300 0.000637

#signif(regionEffect, 3)
print_effect(regionEffect)

#                  slope      se  df       t    pvalue
# grey          0.003560 0.00209 391   1.700  8.94e-02
# turquoise     0.059700 0.00131 399  45.400 1.22e-159
# blue         -0.060900 0.00156 391 -39.000 9.46e-137
# brown         0.004920 0.00292 401   1.680  9.30e-02
# yellow       -0.029500 0.00261 377 -11.300  1.01e-25
# green         0.026600 0.00260 366  10.200  1.10e-21
# red          -0.010800 0.00230 380  -4.680  4.03e-06
# black        -0.000501 0.00284 341  -0.176  8.60e-01
# pink          0.012900 0.00286 382   4.530  8.01e-06
# magenta       0.015300 0.00264 409   5.800  1.35e-08
# purple        0.018900 0.00285 408   6.630  1.05e-10
# greenyellow  -0.018700 0.00301 405  -6.210  1.32e-09
# tan           0.010600 0.00273 371   3.900  1.16e-04
# salmon       -0.027600 0.00186 392 -14.900  5.31e-40
# cyan         -0.036700 0.00262 412 -14.000  1.23e-36
# midnightblue  0.012700 0.00269 402   4.730  3.08e-06
# lightcyan    -0.028900 0.00288 405 -10.100  2.12e-21
# grey60        0.022800 0.00281 391   8.120  6.16e-15


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


print_effect(intEffect)



##################
# make boxplots ##
lab = paste0(substr(rse_gene$BrainRegion,1,4), ":", rse_gene$Dx)
table(lab)
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
# > colnames(m)  WILL CHANGE 
# [1] "(Intercept)"               "DxControl"
# [3] "BrainRegionsACC"           "DxControl:BrainRegionsACC"


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

magmaSig = magma[which(magma$P_JOINT
## line up
mm = match(rowData(rse_gene)$Symbol, magma$SYMBOL)
magGenes = split(magma$SYMBOL[mm], net$colorsLab)

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
