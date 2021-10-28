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
library(dplyr)
library(gprofiler2)


#rm(list = ls())

capabilities()

## load data

load(here('exprs_cutoff', 'rse_gene.Rdata'), verbose = TRUE)
load(here('data', 'degradation_rse_MDDseq_BiPSeq_BothRegions.Rdata'), verbose = TRUE)

table(colData(rse_gene)$PrimaryDx, colData(rse_gene)$BrainRegion)

table(colData(cov_rse)$PrimaryDx, colData(cov_rse)$BrainRegion)
## checks
identical(colnames(rse_gene), colnames(cov_rse)) # TRUE

rse_gene <- rse_gene[, rse_gene$PrimaryDx %in% c("Control", "MDD")]
cov_rse <- cov_rse[, cov_rse$PrimaryDx %in% c("Control", "MDD")]

rse_gene$Dx <- droplevels(rse_gene$PrimaryDx)
cov_rse$Dx <- droplevels(cov_rse$PrimaryDx)

table(colData(rse_gene)$PrimaryDx, colData(rse_gene)$BrainRegion)
table(colData(rse_gene)$Dx, colData(rse_gene)$BrainRegion)


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


table(colData(rse_gene)$PrimaryDx)
table(droplevels(colData(rse_gene)$PrimaryDx))


rse_gene_sACC <- rse_gene[, rse_gene$BrainRegion %in% c("sACC")]
cov_rse_sACC <- cov_rse[, cov_rse$BrainRegion %in% c("sACC")]


table(rse_gene_sACC$BrainRegion)
table(rse_gene_sACC$PrimaryDx)
table(cov_rse_sACC$PrimaryDx)

rse_gene_sACC$PrimaryDx <- droplevels(rse_gene_sACC$PrimaryDx)
cov_rse_sACC$PrimaryDx <- droplevels(cov_rse_sACC$PrimaryDx)

#    MDD Control 
#    228     200 

rse_gene_sACC$PrimaryDx<-relevel(rse_gene_sACC$PrimaryDx, "Control")
cov_rse_sACC$PrimaryDx <- relevel(cov_rse_sACC$PrimaryDx, "Control")


#compute RPKM
assays(rse_gene_sACC)$rpkm = recount::getRPKM(rse_gene_sACC, 'Length')

#######end of model used 

##########
## model #
##########

### put back ERCCsumLogErr term
modJoint_sACC = model.matrix(~PrimaryDx + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 +
                          mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr,
                            data=colData(rse_gene_sACC))


table(modJoint_sACC[, 2])
#  0   1 
#200 228 
## counts from degrafation into log2 scale
degExprs = log2(assays(cov_rse_sACC)$count+1)
##SVA command
k = sva::num.sv(degExprs, modJoint_sACC)
print(k) # 19
### compute principal components from degratation matrix , keep 22 
qSV_mat = prcomp(t(degExprs))$x[,1:k]

colnames(modJoint_sACC)

#  [1] "(Intercept)"       "PrimaryDxMDD"         "AgeDeath"         
#  [4] "SexM"              "snpPC1"            "snpPC2"           
#  [7] "snpPC3"            "mitoRate"          "rRNA_rate"        
# [10] "totalAssignedGene" "RIN"               "ERCCsumLogErr"    

modQsva_sACC = cbind(modJoint_sACC, qSV_mat)

colnames(modQsva_sACC)
# 
#  [1] "(Intercept)"       "PrimaryDxMDD"         "AgeDeath"         
#  [4] "SexM"              "snpPC1"            "snpPC2"           
#  [7] "snpPC3"            "mitoRate"          "rRNA_rate"        
# [10] "totalAssignedGene" "RIN"               "ERCCsumLogErr"    
# [13] "PC1"               "PC2"               "PC3"              
# [16] "PC4"               "PC5"               "PC6"              
# [19] "PC7"               "PC8"               "PC9"              
# [22] "PC10"              "PC11"              "PC12"             
# [25] "PC13"              "PC14"              "PC15"             
# [28] "PC16"              "PC17"              "PC18"             
# [31] "PC19"  


## clean expression
geneExprs = log2(recount::getRPKM(rse_gene_sACC, "Length")+1)

### regress out after variable 5 (protect 1,2,3,4, 5)
geneExprsClean = cleaningY(geneExprs, modQsva_sACC, P=2)

#load(here('wgcna', 'geneExprsClean.rda'), verbose = TRUE)
          
          ##### Load rse data, examine ####

#######end of model used 

## load
load("rdas/constructed_network_signed_bicor.rda", verbose=TRUE)
# Loading objects:
#   net_list
#   net
#   fNames

#net is the main blockwide module command

#net_amyg = net_list[[1]]
net_sACC = net_list[[2]]

net = net_sACC


### added code to create dendogram from WGCNA Tutorial
mergedColors = labels2colors(net$colors)

pdf("dendrogram_102521_sACC.pdf")

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
# [1] 25  4 SACC
### WILL CHANGE WITH AGE PROTECTED
colorDat
#                    num           col Label numGenes
# ENSG00000227232.5    0          grey   ME0    14678
# ENSG00000240731.1    1     turquoise   ME1     1842
# ENSG00000175756.13   2          blue   ME2     1472
# ENSG00000187961.13   3         brown   ME3     1401
# ENSG00000142609.17   4        yellow   ME4      838
# ENSG00000142583.17   5         green   ME5      721
# ENSG00000078369.17   6           red   ME6      641
# ENSG00000162576.16   7         black   ME7      527
# ENSG00000228794.8    8          pink   ME8      461
# ENSG00000162408.10   9       magenta   ME9      434
# ENSG00000188976.10  10        purple  ME10      426
# ENSG00000107404.18  11   greenyellow  ME11      271
# ENSG00000197312.11  12           tan  ME12      236
# ENSG00000158816.15  13        salmon  ME13      223
# ENSG00000121310.16  14          cyan  ME14      214
# ENSG00000116151.13  15  midnightblue  ME15      205
# ENSG00000142634.12  16     lightcyan  ME16      121
# ENSG00000117154.11  17        grey60  ME17      105
# ENSG00000162493.16  18    lightgreen  ME18       98
# ENSG00000074800.13  19   lightyellow  ME19       70
# ENSG00000180287.16  20     royalblue  ME20       53
# ENSG00000117758.13  21       darkred  ME21       47
# ENSG00000142627.12  22     darkgreen  ME22       46
# ENSG00000116731.22  23 darkturquoise  ME23       46
# ENSG00000160058.18  24      darkgrey  ME24       36

# ######################
## gene ontology
### plit genes by each color 
gList = split(rowData(rse_gene_sACC)$EntrezID, 
	factor(net$colorsLab,levels = colorDat$col))
gList = lapply(gList, function(x) as.character(x[!is.na(x)]))
univ = rowData(rse_gene_sACC)$EntrezID
univ = as.character(univ[!is.na(univ)])

# module_names <- c("black", "blue", "brown", "cyan", "green", "greenyellow", "grey", "lightcyan","magenta", "midnightblue", "pink", "purple", "red", "salmon", "tan", "turquoise", "yellow")


#FDR significant
("darkgrey", "green", "turquoise", "magenta", "blue", "pink", "grey60", "royalblue")

#bonf significant
("darkgrey", "green", "turquoise", "magenta")

######################
## gene ontology
### plit genes by each color 
gList = split(rowData(rse_gene_sACC)$EntrezID, 
	factor(net$colorsLab,levels = colorDat$col))
gList = lapply(gList, function(x) as.character(x[!is.na(x)]))
univ = rowData(rse_gene_sACC)$EntrezID
univ = as.character(univ[!is.na(univ)])

go = compareCluster(gList, univ = univ,
	OrgDb = "org.Hs.eg.db", ont = "ALL",
	readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 1)
save(go, file = "rdas/go_enrichment_MDD_Control_sACC_102511_wgcna.rda")
######

goDf = as.data.frame(go)
## below specific for MDD classes 
#goCheck = goDf[goDf$Cluster %in% c("red", "pink", "magenta", "royalblue") &
### changed to include nominally significant categories for the association with depression
goCheck = goDf[goDf$Cluster %in% c("darkgrey", "green", "turquoise", "magenta", "blue", "pink", "grey60", "royalblue") & goDf$qvalue < 0.05,]
goCheck  =goCheck[order(goCheck$pvalue),]
write.csv(goCheck, file = "go_enrichment_MDD_Control_sACC_102511_wgcna_CandidateModules.csv")
##############


##############
## associate eigengenes with brain region
m = modQsva_sACC[,1:2] # this is what was protected
dim(m)
colnames(m)
#### look up MEs in WGCNA 
MEs = net$MEs
## same order 
colnames(MEs) = colorDat$col[match(colnames(MEs), colorDat$Label)]
MEs = MEs[,colorDat$col] #change order of columns
dim(MEs)
#428 25
## check## lmer intercept -- random effect -- check what -1 means 

unique(length(rse_gene_sACC$BrNum))

#statList = lapply(MEs, function(x) summary(lmer(x ~ m + (1|rse_gene_sACC$BrNum) - 1))$coef)

statList = lapply(MEs, function(x) summary(lm(x ~ m))$coeff)

MDDEffect = as.data.frame(t(sapply(statList, function(x) x[2,])))
    
dim(MDDEffect)

colnames(MDDEffect)= c("slope", "se", "t", "pvalue")

print_effect <- function(x) { signif(cbind(x, FDR = p.adjust(x$pvalue, 'fdr')), 3) }
print_effect(MDDEffect)

#                   slope      se       t   pvalue      FDR
# grey           0.001810 0.00469  0.3860 7.00e-01 7.95e-01
# turquoise     -0.018300 0.00461 -3.9600 8.72e-05 7.27e-04
# blue           0.013600 0.00465  2.9200 3.70e-03 1.85e-02
# brown         -0.010300 0.00467 -2.2100 2.80e-02 6.44e-02
# yellow         0.010100 0.00467  2.1700 3.09e-02 6.44e-02
# green         -0.026400 0.00452 -5.8400 1.02e-08 1.27e-07
# red            0.007110 0.00468  1.5200 1.29e-01 2.02e-01
# black         -0.004910 0.00469 -1.0500 2.95e-01 3.89e-01
# pink           0.012200 0.00466  2.6200 9.18e-03 3.28e-02
# magenta       -0.014800 0.00464 -3.1800 1.56e-03 9.78e-03
# purple         0.001260 0.00469  0.2680 7.89e-01 8.22e-01
# greenyellow    0.010200 0.00467  2.1800 3.01e-02 6.44e-02
# tan           -0.001470 0.00469 -0.3130 7.54e-01 8.20e-01
# salmon        -0.005680 0.00469 -1.2100 2.26e-01 3.14e-01
# cyan           0.006120 0.00468  1.3100 1.92e-01 2.83e-01
# midnightblue   0.009420 0.00467  2.0200 4.44e-02 7.40e-02
# lightcyan     -0.002490 0.00469 -0.5310 5.96e-01 7.09e-01
# grey60         0.012200 0.00466  2.6200 9.15e-03 3.28e-02
# lightgreen     0.009460 0.00467  2.0300 4.34e-02 7.40e-02
# lightyellow    0.010500 0.00467  2.2600 2.45e-02 6.44e-02
# royalblue     -0.011900 0.00466 -2.5500 1.12e-02 3.51e-02
# darkred        0.004300 0.00469  0.9170 3.59e-01 4.49e-01
# darkgreen      0.009820 0.00467  2.1000 3.60e-02 6.92e-02
# darkturquoise  0.000124 0.00469  0.0265 9.79e-01 9.79e-01
# darkgrey       0.040400 0.00427  9.4800 1.75e-19 4.38e-18








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
#unclear why you need P values and only significant genes are available 
#magma$FDR_JOINT = p.adjust(magma$P_JOINT, "fdr")
#magma$BONF_JOINT = p.adjust(magma$P_JOINT, "bonf")
magma$Module = net$colorsLab[match(magma$GENE, rowData(rse_gene_sACC)$EntrezID)]
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

magmaSig = magma[magma$BONF_JOINT < 0.05 & magma$GENE %in% rowData(rse_gene_sACC)$EntrezID,]

entrezIDs = split(rowData(rse_gene_sACC)$EntrezID, net$colorsLab)
entrezOverlap = sapply(entrezIDs, function(x) sum(x %in% magmaSig$GENE,na.rm=TRUE))
entrezPresent = sapply(entrezIDs, function(x) sum(x %in% magma$GENE,na.rm=TRUE))
x = data.frame(numOverlap = entrezOverlap, numPresent = entrezPresent)
x = x[colorDat$col,]
x$numGenes = colorDat$numGenes
x$totGene = nrow(rse_gene_sACC)

## pink enrich
mat = matrix(c(5,198, 128, 13894 -5 - 198-128), nr=2)
chisq.test(mat)
getOR(mat)

magmaSig = magma[which(magma$P_JOINT
## line up
mm = match(rowData(rse_gene_sACC)$Symbol, magma$SYMBOL)
magGenes = split(magma$SYMBOL[mm], net$colorsLab)



#####extract genes in module for visualization 

TOM <- get(load("rdas/wgcna_signed_TOM-block.1.RData"))
TOM.mat = as.matrix(TOM)

#TOM = TOMsimilarityFromExpr(t(geneExprsClean), power = 10)

dim(TOM.mat)    
# choose module
module = "lightcyan"
# get list of genes

load("/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/wgcna/geneExprsClean.rda", verbose = TRUE)
dim(geneExprsClean)
probes = rownames(geneExprsClean)

length(probes)
#25212


inModule = (net$colorsLab==module)

modProbes = probes[inModule]
# Select the corresponding Topological Overlap
modTOM = TOM.mat[inModule, inModule]

kIN=softConnectivity(t(geneExprsClean))[, ]

Error in softConnectivity(t(geneExprsClean))[, modProbes] : 
  incorrect number of dimensions
> 

selectHubs = (rank (-kIN) 
vis = exportNetworkToVisANT(modTOM[selectHubs,selectHubs],
file=paste("VisANTInput-", module, "-top30.txt", sep=""),
weighted=TRUE, threshold = 0, probeToGene=
data.frame((GeneAnnotation)$substanceBXH,GeneAnnotation$gene_symbol))

names(colData(rse_gene_sACC))

    
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
file = paste("VisANTInput-", module, ".txt", sep=""),
weighted = TRUE,
threshold = 0)


probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )

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
