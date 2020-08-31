###
library(SummarizedExperiment)
library(recount)
library(edgeR)
library(limma)
library(jaffelab)
library(RColorBrewer)
library(rtracklayer)
library(readxl)
library(clusterProfiler)
library(org.Mm.eg.db)

dir.create("tables")

### load data
load("../data/rse_gene_raw_GoesZandi_n1140.rda",verbose = TRUE)

## make group for plots
rse_gene$Label = paste0(rse_gene$BrainRegion, "_", rse_gene$PrimaryDx)
rse_gene$Label = factor(rse_gene$Label)

# ### make into mouse + h_21 object
# rse_gene = rse_gene[rowData(rse_gene)$gene_species == "mouse" | 
# 	seqnames(rse_gene) == "human_chr21",]

## metric boxplot to drop sample
pdf("pdfs/metrics_vs_grouplabels_predrop.pdf")
par(mar=c(13,6,2,2), cex.axis=1.8,cex.lab=1.8)
palette(brewer.pal(4,"Set1"))
boxplot(mitoRate ~ Label, data = colData(rse_gene), 
        xlab="",outline = FALSE, ylim = range(rse_gene$mitoRate),
        ylab = "chrM mapping rate", las=3)
text(jitter(as.numeric(rse_gene$Label),amount=0.15),
     rse_gene$mitoRate,letters[as.numeric(factor(rse_gene$BrNum))], 
     col=as.numeric(rse_gene$BrainRegion),font=2)
abline(h=0.5, lty=2)
boxplot(overallMapRate ~ Label, data = colData(rse_gene), 
        xlab="",outline = FALSE, ylim = range(rse_gene$overallMapRate),
        ylab = "overall mapping rate", las=3)
text(jitter(as.numeric(rse_gene$Label),amount=0.15),
     rse_gene$overallMapRate,letters[as.numeric(factor(rse_gene$BrNum))], 
     col=as.numeric(rse_gene$BrainRegion),font=2)
abline(h=0.86, lty=2)
dev.off()

## filter for expressed
gene_rpkm = getRPKM(rse_gene,"Length")
rse_gene_filter = rse_gene[rowMeans(gene_rpkm) > 0.1,]
gene_rpkm_filter = gene_rpkm[rowMeans(gene_rpkm) > 0.1,]
table(droplevels(seqnames(rse_gene_filter)))

## drop one sample
keepIndex=rse_gene_filter$overallMapRate > 0.86
gene_rpkm_filter = gene_rpkm_filter[,keepIndex]
rse_gene_filter = rse_gene_filter[,keepIndex]
dim(rse_gene_filter)

## get expression
geneExprs_filter = log2(gene_rpkm_filter+1)

#####################
#### explore pcs
pca = prcomp(t(geneExprs_filter))
pca_vars = getPcaVars(pca)
pca_vars_lab = paste0("PC", seq(along=pca_vars), ": ",
                      pca_vars, "% Var Expl")

#### pc1v2
pdf("pdfs/pc1_vs_2.pdf")
par(mar=c(5,6,2,2),cex.axis=1.8,cex.lab=1.8)
palette(brewer.pal(6,"Set1"))
plot(pca$x, pch=as.numeric(rse_gene_filter$Label) + 20, 
     bg=rse_gene_filter$Label,cex=1.5,
     xlab = pca_vars_lab[1], ylab = pca_vars_lab[2])
legend("topright", levels(rse_gene_filter$Label),
       col = 1:4,	pch = 19, cex=1.5)
# legend("topright", levels(rse_gene_filter$BrainRegion),
# 	col = 1:4,	pch = 19, cex=1.5)
# legend("bottomright", levels(rse_gene_filter$Group),
# 	pt.bg = "black",pch = 21:22,cex=1.5)
dev.off()

#### PCs vs stuff

pdf("pdfs/pcs_vs_grouplabels.pdf")
par(mar=c(13,6,2,2), cex.axis=1.8,cex.lab=1.8)
palette(brewer.pal(4,"Set1"))
for(i in 1:10) {
  boxplot(pca$x[,i] ~ rse_gene_filter$Label, xlab="",
          outline = FALSE, ylim = range(pca$x[,i]),
          ylab = pca_vars_lab[i], las=3)
  text(jitter(as.numeric(rse_gene_filter$Label),amount=0.15),pca$x[,i],
       letters[as.numeric(factor(rse_gene_filter$BrNum))], 
       col=as.numeric(rse_gene_filter$BrainRegion),font=2)
}
dev.off()

### variables
rse_gene_filter$BrNum = factor(rse_gene_filter$BrNum)

vars = c("rna_conc_ug_u_l", "x260_280", "total_rna_given_ug",
         "overallMapRate", "overallMapRate", 
         "mitoRate", "totalAssignedGene",
         "totalUnassignedAmbig", "totalUnassignedMulti",	"rRNA_rate")
names(vars) = c("RNA Conc", "260/280", "Total RNA", "Overall Map Rate",
                "Corcordant Map Rate", "chrM Map Rate", "Exonic Assignment Rate",
                "Ambig Assign Rate", "Multi Assign Rate", "rRNA Assignment Rate")

## correlation
ccPcaMetrics = cor(pca$x[,1:20], 
                   as.data.frame(colData(rse_gene_filter)[,vars]))
tPcaMetrics = ccPcaMetrics/sqrt((1-ccPcaMetrics^2)/(ncol(rse_gene_filter)-2))
pvalPcaMetrics = 2*pt(-abs(tPcaMetrics), df = ncol(rse_gene_filter)-1)

## add anova metrics
regVars = c("BrainRegion", "Group","BrNum")
names(regVars) = c("BrainRegion", "MAC21", "Mouse")
anovaPval = sapply(regVars, function(x) {
  y = colData(rse_gene_filter)[,x]
  apply(pca$x[,1:20], 2, function(p) {
    anova(lm(p ~ y))$`Pr(>F)`[1]
  })
})

library(lattice)
pdf("pdfs/pcs_vs_metrics_heatmap.pdf",w=9)
colnames(pvalPcaMetrics) = names(vars)
logPvalMat = -log10(cbind(pvalPcaMetrics, anovaPval))
logPvalMat[logPvalMat>16] = 16
theSeq = seq(0,16,by=0.1)
my.col <- colorRampPalette(c("white","blue"))(length(theSeq))
print(levelplot(logPvalMat, aspect = "fill", 
                at = theSeq,pretty=TRUE,xlab="",ylab="",
                scales=list(x=list(rot=90, cex=1.2), y=list(cex=1.2)),
                panel = panel.levelplot.raster, col.regions = my.col))
dev.off()

#################
## de modeling ##
#################

## norm counts
dge = DGEList(counts = assays(rse_gene_filter)$counts, 
              genes = rowData(rse_gene_filter))
dge = calcNormFactors(dge)

## split by BrainRegion
tIndexes = splitit(rse_gene_filter$BrainRegion)

#########
## human modeling first
human_chr21_index = which(seqnames(rse_gene_filter) == "human_chr21")

statList_human_chr21 = lapply(tIndexes, function(ii) {
  ## filter to human
  d = dge[human_chr21_index,ii]
  
  ## mean-variance
  mod = model.matrix(~Group + totalAssignedGene + mitoRate,
                     data=colData(rse_gene_filter[,ii]))
  
  vGene = voom(d,mod,plot=FALSE)
  fitGene= lmFit(vGene)
  ebGene = eBayes(fitGene)
  tt = topTable(ebGene,coef=2,
                p.value = 1,number=nrow(d),sort="none")
  tt[,c( "logFC", "AveExpr", "t", "P.Value", "adj.P.Val","B")]
})
stats_human_chr21 = do.call("cbind", statList_human_chr21)
colMeans(stats_human_chr21[,grep("P.Val", colnames(stats_human_chr21))] < 0.05)

### mouse next
mouse_chrs_index  = which(seqnames(rse_gene_filter) != "human_chr21")


statList_mouse = lapply(tIndexes, function(ii) {
  ## filter to human
  d = dge[mouse_chrs_index,ii]
  
  ## mean-variance
  mod = model.matrix(~Group + totalAssignedGene + mitoRate,
                     data=colData(rse_gene_filter[,ii]))
  
  vGene = voom(d,mod,plot=FALSE)
  fitGene= lmFit(vGene)
  ebGene = eBayes(fitGene)
  tt = topTable(ebGene,coef=2,
                p.value = 1,number=nrow(d),sort="none")
  tt[,c( "logFC", "AveExpr", "t", "P.Value", "adj.P.Val","B")]
})
stats_mouse = do.call("cbind", statList_mouse)
colMeans(stats_mouse[,grep("P.Val", colnames(stats_mouse))] < 0.05)
colSums(stats_mouse[,grep("adj", colnames(stats_mouse))] < 0.1)

## one overall model for mouse
mod = model.matrix(~Group*BrainRegion + totalAssignedGene + mitoRate,
                   data=colData(rse_gene_filter))
vGene = voom(dge[mouse_chrs_index,],mod,plot=FALSE)
fitGene= lmFit(vGene)
ebGene = eBayes(fitGene)

outGene_MAC21 = topTable(ebGene,coef=2,
                         p.value = 1,number=nrow(dge[mouse_chrs_index,]),sort="none")
table(	outGene_MAC21$adj.P.Val < 0.05)

outGene_MAC21_fatEffect = topTable(ebGene,coef=8:10,
                                   p.value = 1,number=nrow(dge[mouse_chrs_index,]),sort="none")
table(outGene_MAC21_fatEffect$adj.P.Val < 0.05,
      outGene_MAC21$adj.P.Val < 0.05,
      dnn = c("FATxMAC21", "MAC21"))

outGene_MAC21_list = lapply(8:10, function(i) {
  topTable(ebGene,coef=i,	p.value = 1,
           number=nrow(dge[mouse_chrs_index,]),
           sort="none")
})

outGene_MAC21_inter = do.call("cbind", outGene_MAC21_list)