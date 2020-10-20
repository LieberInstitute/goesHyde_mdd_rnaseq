### Combined-batch iPSC>neuron data from Hiler/Maher group, for reference =========
library(RColorBrewer)
library(ggrepel)
pal = brewer.pal(10, "Set3")
pal[2] = "#8B795E" # darken NPCs
load("/dcl01/lieber/ajaffe/Matt/BradyMaher/hiler_Mar-BlueCombined/rdas/rse_gene_n133-Mar-BlueCombined_rpkmFiltered_MTN02Apr2019.rda",
     verbose=T)
# rse_comb
## FIRST remove Blue batch N11C2's - decided with DaHi 15Apr2019
rse_comb <- rse_comb[ ,-which(rse_comb$Line=="N11C2" & rse_comb$ExpBatch=="Blue")]
### cell composition
load("/dcl02/lieber/ajaffe/libd_stem_timecourse/deconvolution/cell_type/singleCell_iPSC_quake_coefEsts_calibration_Zscale.rda",
     verbose=T)
# coefEsts, mergeMarkerMeanExprsZ
geneExprs = log2(assays(rse_comb)$rpkm+1)
table(rownames(coefEsts) %in% rownames(geneExprs))
#FALSE  TRUE
#48   180
geneExprs_scale = scale(geneExprs[rownames(coefEsts)[rownames(coefEsts) %in% rownames(geneExprs)], ])
cellPropEsts = minfi:::projectCellType(geneExprs_scale,coefEsts)
colData(rse_comb) = cbind(colData(rse_comb), as.data.frame(cellPropEsts))
## make plots
#pdf("pdfs/estRnaFractions_n127-Mar-BlueCombined_MTN23Apr2019.pdf",h=6,w=8)
par(mar=c(7,6,3,2), cex.axis=0.8,cex.lab=1.5, cex.main=1.5)
for(i in 67:76) {
  boxplot(colData(rse_comb)[ ,i] ~  as.numeric(ss(rse_comb$DiffDay, "Day", 2)) +
            relevel(factor(rse_comb$Genotype),"WT"), las=3,
          outline = FALSE, ylab = "RNA Fraction",
          main = paste0(colnames(colData(rse_comb))[i], " Composition by Fraction")
  )
  mtext(text="Differentiation Day & Genotype",side=1, padj=6.5, font=2)
}
#dev.off()
cellPropEsts_scaled <- prop.table(cellPropEsts,1)
gIndexes = splitit(factor(rse_comb$DiffDay,levels=c("Day0","Day3","Day7","Day10","Day17","Day20")))
cellPropEsts_groupMeans = sapply(gIndexes, function(ii) colMeans(cellPropEsts[ii,]))
cellPropEsts_groupSEs = sapply(gIndexes, function(ii) apply(cellPropEsts[ii,],2,function(x) sd(x)/sqrt(length(x))))
## Line plot - didn't do MTN 25Feb2019
#pdf("pdfs/cellTypeDecon_lineplot_Hiler_n63_MTN04Oct2018.pdf",h=5,w=5)
palette(pal)
par(mar=c(10,6,2,2),cex.axis=1.2,cex.lab=1.3)
plot(cellPropEsts_groupMeans[1,], type="b",xaxt="n",
     lty=1, pch=19,lwd=2, ylim = c(0,1),col=1,xlab="",
     ylab = "Class Proportion",cex=1.3)
axis(1,at=seq(along=gIndexes),gsub(" ", "", colnames(cellPropEsts_groupMeans)), 
     las=3,cex.axis=1.5)
segments(x0=seq(along=gIndexes), x1=seq(along=gIndexes), col=1,lwd=2, 
         y0=cellPropEsts_groupMeans[1,] - 2*cellPropEsts_groupSEs[1,], 
         y1=cellPropEsts_groupMeans[1,] + 2*cellPropEsts_groupSEs[1,])
for(i in 2:10) {
  lines(cellPropEsts_groupMeans[i,], type="b",
        lty=1, pch=19,lwd=2, col=i,cex=1.3)
  segments(x0=seq(along=gIndexes), x1=seq(along=gIndexes), col=i,lwd=2, 
           y0=cellPropEsts_groupMeans[i,] - 2*cellPropEsts_groupSEs[i,], 
           y1=cellPropEsts_groupMeans[i,] + 2*cellPropEsts_groupSEs[i,])
}
#dev.off()	
## Bar plot
#pdf("pdfs/cellTypeDecon_barplot_Hiler_n127-Mar-BlueCombined_MTN23Apr2019.pdf",w=20,h=8)
palette(pal)
par(mar = c(14,6,2,2), cex.axis=1.5,cex.lab=2)
oo= order(factor(rse_comb$DiffDay,levels=c("Day0","Day3","Day7","Day10","Day17","Day20")),
          relevel(factor(rse_comb$Genotype), "WT"))
bp = barplot(t(cellPropEsts_scaled[oo,]), col = pal,
             ylim = c(0,1.45),xaxt = "n", yaxt= "n",ylab="Class Proportion            ")
g = split(bp, factor(rse_comb$DiffDay,levels=c("Day0","Day3","Day7","Day10","Day17","Day20"))[oo])
axis(1,at=sapply(g,mean), font=2, labels = names(g),las=3, cex.axis=2)
abline(v = sapply(g,min)[-1]-0.6,lwd=3.5,col="black")
legend("top", colnames(cellPropEsts_scaled), pch = 15, 
       col = 1:10,bg="white", nc = 5, cex=2, pt.cex=3)
axis(2, at=seq(0,1,by=0.25))
# Annotate with Genotype identifiers
text(x=unname(unlist(g)),y=1.03, labels=ifelse(rse_comb$Genotype[oo]=="WT","W","P"),
     cex=0.9, font=2, col=c(ifelse(rse_comb$Genotype[oo]=="WT","black","red")))
#dev.off()
## END REFERENCE =========
