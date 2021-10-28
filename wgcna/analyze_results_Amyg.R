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

capabilities()

#rm(list = ls())

## load data

load(here('exprs_cutoff', 'rse_gene.Rdata'), verbose = TRUE)
load(here('data', 'degradation_rse_MDDseq_BiPSeq_BothRegions.Rdata'), verbose = TRUE)

table(colData(rse_gene)$PrimaryDx, colData(rse_gene)$BrainRegion)

table(colData(cov_rse)$PrimaryDx, colData(cov_rse)$BrainRegion)
## checks
identical(colnames(rse_gene), colnames(cov_rse)) # TRUE

rse_gene <- rse_gene[, rse_gene$PrimaryDx %in% c("Control", "MDD")]
cov_rse <- cov_rse[, cov_rse$PrimaryDx %in% c("Control", "MDD")]

rse_gene$PrimaryDx <- droplevels(rse_gene$PrimaryDx)
cov_rse$PrimaryDx <- droplevels(cov_rse$PrimaryDx)

table(colData(rse_gene)$PrimaryDx, colData(rse_gene)$BrainRegion)
table(colData(rse_gene)$PrimaryDx, colData(rse_gene)$BrainRegion)

###########


table(colData(rse_gene)$PrimaryDx)
table(droplevels(colData(rse_gene)$PrimaryDx))


#rse_gene_sACC <- rse_gene[, rse_gene$BrainRegion %in% c("sACC")]
#cov_rse_sACC <- cov_rse[, cov_rse$BrainRegion %in% c("sACC")]


rse_gene_Amyg <- rse_gene[, rse_gene$BrainRegion %in% c("Amygdala")]
cov_rse_Amyg <- cov_rse[, cov_rse$BrainRegion %in% c("Amygdala")]


table(colData(rse_gene)$PrimaryDx)
table(droplevels(colData(rse_gene)$PrimaryDx))


rse_gene_Amyg <- rse_gene[, rse_gene$BrainRegion %in% c("Amygdala")]
cov_rse_Amyg <- cov_rse[, cov_rse$BrainRegion %in% c("Amygdala")]


table(rse_gene_Amyg$BrainRegion)
table(rse_gene_Amyg$PrimaryDx)
table(cov_rse_Amyg$PrimaryDx)

rse_gene_Amyg$PrimaryDx <- droplevels(rse_gene_Amyg$PrimaryDx)
cov_rse_Amyg$PrimaryDx <- droplevels(cov_rse_Amyg$PrimaryDx)

rse_gene_Amyg$PrimaryDx<-relevel(rse_gene_Amyg$PrimaryDx, "Control")
cov_rse_Amyg$PrimaryDx <- relevel(cov_rse_Amyg$PrimaryDx, "Control")

#Control     MDD 
#    187     231 
 

#compute RPKM
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')

#######end of model used 

##########
## model #
##########

### put back ERCCsumLogErr term
modJoint_Amyg = model.matrix(~PrimaryDx + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 +
                          mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr,
                            data=colData(rse_gene_Amyg))

table(modJoint_Amyg[, 2])


## counts from degrafation into log2 scale
degExprs = log2(assays(cov_rse_Amyg)$count+1)
##SVA command
k = sva::num.sv(degExprs, modJoint_Amyg)
print(k) # 18
### compute principal components from degratation matrix , keep 22 
qSV_mat = prcomp(t(degExprs))$x[,1:k]

colnames(modJoint_Amyg)

#  [1] "(Intercept)"       "PrimaryDxMDD"         "AgeDeath"         
#  [4] "SexM"              "snpPC1"            "snpPC2"           
#  [7] "snpPC3"            "mitoRate"          "rRNA_rate"        
# [10] "totalAssignedGene" "RIN"               "ERCCsumLogErr"    

modQsva_Amyg = cbind(modJoint_Amyg, qSV_mat)

colnames(modQsva_Amyg)
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
# # > 

## clean expression
geneExprs = log2(recount::getRPKM(rse_gene_Amyg, "Length")+1)

### regress out after variable 5 (protect 1,2,3,4, 5)
geneExprsClean = cleaningY(geneExprs, modQsva_Amyg, P=2)

#load(here('wgcna', 'geneExprsClean.rda'), verbose = TRUE)
          
          ##### Load rse data, examine ####

#######end of model used 

## load
load("wgcna/rdas/constructed_network_signed_bicor.rda", verbose=TRUE)
# Loading objects:
#   net_list
#   net
#   fNames

#net is the main blockwide module command

net_amyg = net_list[[1]]

net = net_amyg


### added code to create dendogram from WGCNA Tutorial
mergedColors = labels2colors(net$colors)

pdf("dendrogram_102321_Amyg.pdf")

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
# [1] 33  4 Amyg
### WILL CHANGE WITH AGE PROTECTED
colorDat

#                    num           col Label numGenes
# ENSG00000227232.5    0          grey   ME0    15872
# ENSG00000240731.1    1     turquoise   ME1     1196
# ENSG00000142583.17   2          blue   ME2      838
# ENSG00000142609.17   3         brown   ME3      666
# ENSG00000188157.13   4        yellow   ME4      606
# ENSG00000197921.5    5         green   ME5      564
# ENSG00000236948.2    6           red   ME6      533
# ENSG00000011021.21   7         black   ME7      473
# ENSG00000196581.10   8          pink   ME8      458
# ENSG00000228794.8    9       magenta   ME9      435
# ENSG00000107404.18  10        purple  ME10      404
# ENSG00000230415.1   11   greenyellow  ME11      381
# ENSG00000224051.6   12           tan  ME12      306
# ENSG00000217801.9   13        salmon  ME13      302
# ENSG00000188976.10  14          cyan  ME14      290
# ENSG00000187608.8   15  midnightblue  ME15      273
# ENSG00000171824.13  16     lightcyan  ME16      219
# ENSG00000187634.11  17        grey60  ME17      165
# ENSG00000162591.15  18    lightgreen  ME18      146
# ENSG00000157933.9   19   lightyellow  ME19      129
# ENSG00000178965.13  20     royalblue  ME20      115
# ENSG00000116329.10  21       darkred  ME21      114
# ENSG00000173406.15  22     darkgreen  ME22      100
# ENSG00000189337.16  23 darkturquoise  ME23       99
# ENSG00000186094.16  24      darkgrey  ME24       93
# ENSG00000183114.7   25        orange  ME25       69
# ENSG00000162493.16  26    darkorange  ME26       67
# ENSG00000183888.4   27         white  ME27       57
# ENSG00000184007.18  28       skyblue  ME28       54
# ENSG00000169641.13  29   saddlebrown  ME29       52
# ENSG00000171502.14  30     steelblue  ME30       47
# ENSG00000142627.12  31 paleturquoise  ME31       46
# ENSG00000070886.11  32        violet  ME32       43
> 
# ######################


geneinfo <- as.data.frame(rowRanges(rse_gene_Amyg))
str(geneinfo)
geneinfo$Module = net$colorsLab
table(geneinfo$Module)

write.csv(geneinfo, file = "geneinfo_Amyg.csv")
MEs = net$MEs
## same order 
colnames(MEs) = colorDat$col[match(colnames(MEs), colorDat$Label)]
MEs = MEs[,colorDat$col] #change order of columns
pheno <- data.frame(colData(rse_gene_Amyg), MEs)
write.csv(pheno, file = "colData_ME_Amyg.csv")



#### split genes with enselble ID 
gList = split(rowData(rse_gene_Amyg)$ensemblID, 
	factor(net$colorsLab,levels = colorDat$col))
gList = lapply(gList, function(x) as.character(x[!is.na(x)]))
univ = rowData(rse_gene_Amyg)$ensemblID
univ = as.character(univ[!is.na(univ)])


gostres_multi <-gost(
query = gList,
organism = "hsapiens",
ordered_query = FALSE,
multi_query = FALSE,
significant = TRUE,
exclude_iea = FALSE,
measure_underrepresentation = FALSE,
evcodes = TRUE,
user_threshold = 0.05,
correction_method = "fdr",
domain_scope = c("custom"),
custom_bg = univ,
numeric_ns = "",
sources = "GO",
as_short_link = FALSE
)

save(gostres_multi, file = "gprofiler_results_Amygdala_102621.Rda")


dim(gostres_multi$result)
#head(gostres_multi$result)
names(gostres_multi$result)
col_interest <- c("query","source", "term_name","term_size","query_size","intersection_size","p_value"  )
head(gostres_multi$result[0:12][, col_interest])  

df <- (gostres_multi$result[0:12][, col_interest])  
df$logFDR <- -log10(df$p_value) 
glimpse(df)
df <- as_tibble(df)

write.table(gostres_multi$result[0:13], file = "multi_query_resuls_amygdala_102321.csv", sep =",")

######## 
######################
## gene ontology
gList = split(rowData(rse_gene_Amyg)$EntrezID, 
	factor(net$colorsLab,levels = colorDat$col))
gList = lapply(gList, function(x) as.character(x[!is.na(x)]))
univ = rowData(rse_gene_Amyg)$EntrezID
univ = as.character(univ[!is.na(univ)])

#print list of genes in each module
lapply(gList, function(x) write.table( data.frame(x), 'test.csv'  , append= T, sep=',' ))

lapply(1:length(gList), function(x) write.table(t(as.data.frame(gList[x])), 
                                              'test3.csv', append= T, sep=',', 
                                              quote = F, col.names = F))
###


go = compareCluster(gList, univ = univ,
	OrgDb = "org.Hs.eg.db", ont = "ALL",
	readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 1)


save(go, file = "go_enrichment_MDDControl_wgcna_amyg_1023.21.rda")

go <- get(load("go_enrichment_MDDControl_wgcna_amyg_1023.21.rda"))

goDf = as.data.frame(go)
goCheck = goDf[goDf$Cluster %in% c("blue", "blck", "magenta", "lightcyan", "lightgreen", "lightyellow", "darkturquoise") &
	goDf$qvalue < 0.05,]
goCheck  =goCheck[order(goCheck$pvalue),]
write.csv(goCheck, file = "go_enrichment_MDDControl_wgcna_Amyg_CandidateModules.csv")

##############












#####----
head(gost_cyan$result[0:12][, c(9,10,11,4,6,8,7,3)])
df <- gost_cyan$result[0:12][, c(9,10,11,4,6,8,7,3)]
df$logFDR <- -log10(df$p_value) 

gostplot(gostres_multi, capped = FALSE, interactive = TRUE)

p <- gostplot(gost_cyan, capped = FALSE, interactive = FALSE)
pp <- publish_gostplot(p, highlight_terms = c(sig_cat$term_id), 
                       width = NA, height = NA, filename = NULL )
pp

write.table(gostres_multi$result[0:13], file = "multi_query_resuls.csv", sep =",")

gostres$result[1:20, c(1,2,3,4,6,9,10,11)]
###----


## PLOT RESULTS

module_names <- c("black", "blue", "brown", "cyan", "green", "greenyellow", "grey", "lightcyan","magenta", "midnightblue", "pink", "purple", "red", "salmon", "tan", "turquoise", "yellow")
for (i in module_names) {
df_subset <-  df %>% filter(query == i) %>% mutate(term_name = fct_reorder(term_name, logFDR))
        p <- ggplot(df_subset, mapping = aes(x = logFDR, y = term_name)) + 
        geom_bar(stat = "identity", color = "blue3", fill = "blue3", width = 0.5) +
    	labs(title="GO Enrichment",x = "-log10 FDR") +
    	theme(axis.text.x = element_text(angle=0,hjust=1,vjust=1.0,size = 10),
    		axis.text.y = element_text(size = 12),
    		plot.title = element_text(hjust = 0.5)) +
        theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        geom_vline(xintercept=1.3, color = "red", size = 1.5) +
#        scale_x_continuous(limits = c(0,8))
     pdf(paste0(i, "GO_enrichement.pdf"))
    print(p)
dev.off()
    }
    

##############
##############
## associate eigengenes with MDD
m = modQsva_Amyg[,1:2] # this is what was protected
dim(m)
colnames(m)
#### look up MEs in WGCNA 
MEs = net$MEs
## same order 
colnames(MEs) = colorDat$col[match(colnames(MEs), colorDat$Label)]
MEs = MEs[,colorDat$col] #change order of columns
dim(MEs)
#418 33
## check## lmer intercept -- random effect -- check what -1 means 

unique(length(rse_gene_Amyg$BrNum))

#statList = lapply(MEs, function(x) summary(lmer(x ~ m + (1|rse_gene_sACC$BrNum) - 1))$coef)

statList = lapply(MEs, function(x) summary(lm(x ~ m))$coeff)

MDDEffect = as.data.frame(t(sapply(statList, function(x) x[2,])))
    
dim(MDDEffect)

colnames(MDDEffect)= c("slope", "se", "t", "pvalue")

print_effect <- function(x) { signif(cbind(x, FDR = p.adjust(x$pvalue, 'fdr')), 3) }
print_effect(MDDEffect)

statList_Amyg = statList
MDDEffect_Amyg = MDDEffect

#                 slope      se       t   pvalue    FDR
# grey           5.90e-04 0.00482  0.1220 0.903000 0.9510
# turquoise      3.25e-03 0.00482  0.6740 0.500000 0.7760
# blue          -1.32e-02 0.00478 -2.7600 0.006110 0.0301
# brown         -5.64e-03 0.00482 -1.1700 0.242000 0.4660
# yellow         7.99e-03 0.00481  1.6600 0.097200 0.2460
# green          3.17e-03 0.00482  0.6580 0.511000 0.7760
# red           -4.71e-04 0.00482 -0.0976 0.922000 0.9510
# black          1.45e-02 0.00477  3.0400 0.002520 0.0208
# pink           2.57e-03 0.00482  0.5340 0.594000 0.7760
# magenta       -1.63e-02 0.00476 -3.4300 0.000658 0.0128
# purple        -1.53e-03 0.00482 -0.3160 0.752000 0.9150
# greenyellow   -5.37e-03 0.00482 -1.1200 0.265000 0.4660
# tan           -6.81e-03 0.00481 -1.4200 0.158000 0.3470
# salmon        -3.00e-03 0.00482 -0.6220 0.535000 0.7760
# cyan           1.37e-03 0.00482  0.2850 0.776000 0.9150
# midnightblue   1.20e-02 0.00479  2.5000 0.012900 0.0534
# lightcyan     -1.61e-02 0.00476 -3.3900 0.000775 0.0128
# grey60        -2.78e-03 0.00482 -0.5780 0.564000 0.7760
# lightgreen    -1.34e-02 0.00478 -2.8000 0.005280 0.0301
# lightyellow    1.31e-02 0.00478  2.7400 0.006380 0.0301
# royalblue      7.82e-03 0.00481  1.6300 0.105000 0.2460
# darkred        9.97e-03 0.00480  2.0800 0.038300 0.1410
# darkgreen     -2.45e-03 0.00482 -0.5080 0.612000 0.7760
# darkturquoise  1.55e-02 0.00476  3.2500 0.001250 0.0137
# darkgrey       5.63e-03 0.00482  1.1700 0.243000 0.4660
# orange        -5.34e-03 0.00482 -1.1100 0.268000 0.4660
# darkorange     8.31e-03 0.00481  1.7300 0.084700 0.2330
# white         -9.55e-05 0.00482 -0.0198 0.984000 0.9840
# skyblue        2.74e-03 0.00482  0.5690 0.570000 0.7760
# saddlebrown    9.70e-03 0.00480  2.0200 0.044000 0.1450
# steelblue      5.55e-04 0.00482  0.1150 0.909000 0.9510
# paleturquoise  6.75e-04 0.00482  0.1400 0.889000 0.9510
# violet        -9.14e-03 0.00480 -1.9000 0.057800 0.1730

save(statList_Amyg, MDDEffect_Amyg,statList_sACC, MDDEffect_sACC, file = "WGCNA_module_eigengenes_association_by_region.rda")



##################
# make boxplots ##
lab = paste0(substr(rse_gene_Amyg$BrainRegion,1,4), ":", rse_gene_Amyg$PrimaryDx)
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
	boxplot(clean_ME[,i] ~ rse_gene_Amyg$PrimaryDx, outline = FALSE, xlab="",
		ylim = quantile(unlist(clean_ME),c(0.001,0.999)),main = colnames(MEs)[i],
		ylab = "Module Eigengene (Adj)")
	points(clean_ME[,i] ~ jitter(as.numeric(rse_gene_Amyg$PrimaryDx),amount=0.1), pch=21, bg=lab)
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
magma = read_excel("/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/PGC_BP_Magma_table.xlsx", skip = 2)
#updated 
magma = read.table(file = "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/Howard_MDD_MAGMA_genes.csv", sep = ",")
head(magma)

colnames(magma) <- ("ensemblID")


magma = as.data.frame(magma)
#magma$FDR_JOINT = p.adjust(magma$P_JOINT, "fdr")
#magma$BONF_JOINT = p.adjust(magma$P_JOINT, "bonf")
#magma$Module = net$colorsLab[match(magma$GENE, rowData(rse_gene)$EntrezID)]

magma$Module = net$colorsLab[match(magma$ensemblID, rowData(rse_gene_Amyg)$ensemblID)]

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

entrezIDs = split(rowData(rse_gene_Amyg)$EntrezID, net$colorsLab)
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

names(colData(rse_gene))

    
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
