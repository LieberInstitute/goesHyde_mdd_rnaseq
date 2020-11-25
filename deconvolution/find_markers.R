library(SingleCellExperiment)
# library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
# library(batchelor)
# library(DropletUtils)
library(jaffelab)
library(limma)
library(here)

## Load rse_gene data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

#### sACC data ####
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose = TRUE)

sce.sacc <- sce.sacc[,sce.sacc$cellType != "Ambig.lowNtrxts",]
sce.sacc$cellType <- droplevels(sce.sacc$cellType)
## Add cellType.broad
sce.sacc$cellType.Broad <- ss(as.character(sce.sacc$cellType), "\\.", 1)
sce.sacc$cellType.Broad <- factor(sce.sacc$cellType.Broad)
## Match rownames
rownames(sce.sacc) <- rowData(sce.sacc)$ID
table(rownames(sce.sacc) %in% rownames(rse_gene))
# FALSE  TRUE 
# 15021 18517 

#### Filter data ####
## filter for common between sc and bulk data
common_genes <- rownames(rse_gene)[rownames(rse_gene) %in% rownames(sce.sacc)]
sce.sacc <- sce.sacc[common_genes, ]

# Remove 0 genes across all nuclei
sce.sacc <- sce.sacc[!rowSums(assay(sce.sacc, "counts"))==0, ] 
dim(sce.sacc)
# [1] 17785  7004
## Filter

cellSubtype.idx <- splitit(sce.sacc$cellType.Broad)
medianNon0.idx <- lapply(cellSubtype.idx, function(x){
  apply(as.matrix(assay(sce.sacc, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

sapply(medianNon0.idx, sum)
# Astro Excit Inhib Micro Oligo   OPC 
# 1501  6892  5658   992  1259  2279

## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.sacc), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

markers.sacc.t.1vAll <- list()
for(i in levels(sce.sacc$cellType.Broad)){
  # Make temporary contrast
  message(i)
  sce.sacc$contrast <- ifelse(sce.sacc$cellType.Broad==i, 1, 0)
  # Test cluster vs. all
  markers.sacc.t.1vAll[[i]] <- findMarkers(sce.sacc, groups=sce.sacc$contrast,
                                          assay.type="logcounts", design=mod, test="t",
                                          direction="up", pval.type="all", full.stats=T)
}

temp.1vAll <- list()
for(i in levels(sce.sacc$cellType.Broad)){
  # Make temporary contrast
  message(i)
  sce.sacc$contrast <- ifelse(sce.sacc$cellType.Broad==i, 1, 0)
  # Test cluster vs. all
  temp.1vAll[[i]] <- findMarkers(sce.sacc, groups=sce.sacc$contrast,
                                           assay.type="logcounts", design=mod, test="t",
                                           std.lfc=TRUE,
                                           direction="up", pval.type="all", full.stats=T)
}

# Replace that empty slot with the entry with the actual stats
# What is stats 1 vs. stats 0? 
markers.sacc.t.1vAll <- lapply(markers.sacc.t.1vAll, function(x){ x[[2]] })
# Same for that with std.lfc
temp.1vAll <- lapply(temp.1vAll, function(x){ x[[2]] })


# Now just pull from the 'stats.0' DataFrame column
markers.sacc.t.1vAll <- lapply(markers.sacc.t.1vAll, function(x){ x$stats.0 })
temp.1vAll <- lapply(temp.1vAll, function(x){ x$stats.0 })
m <- markers.sacc.t.1vAll
t <- temp.1vAll
# Re-name std.lfc column and add to the first result
for(i in levels(sce.sacc$cellType.Broad)){
  colnames(temp.1vAll[[i]])[1] <- "std.logFC"
  markers.sacc.t.1vAll[[i]] <- cbind(markers.sacc.t.1vAll[[i]], temp.1vAll[[i]]$std.logFC)
  # Oh the colname is kept weird
  colnames(markers.sacc.t.1vAll[[i]])[4] <- "std.logFC"
  # Then re-organize
  markers.sacc.t.1vAll[[i]] <- markers.sacc.t.1vAll[[i]][ ,c("logFC","std.logFC","log.p.value","log.FDR")]
  markers.sacc.t.1vAll[[i]] <- cbind(markers.sacc.t.1vAll[[i]],
                                    medianNon0.idx[[i]][match(rownames(markers.sacc.t.1vAll[[i]]),
                                                              names(medianNon0.idx[[i]]))])
  colnames(markers.sacc.t.1vAll[[i]])[5] <- "non0median"
}

save(markers.sacc.t.1vAll, file = here("deconvolution","data","markers_broad.sacc.Rdata"))

markerList.t.1vAll <- lapply(markers.sacc.t.1vAll, function(x){
  rownames(x)[x$log.FDR < log10(0.000001) & x$non0median==TRUE]
}
)

lengths(markerList.t.1vAll)
# Astro Excit Inhib Micro Oligo   OPC 
# 832  5162  4186   708   841  1160    

genes.top40.t <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(paste0("plots/expr/sACC_t-sn-level_1vALL_top40markers-REFINED-",i,"_logExprs.png"), height=1900, width=1200)
  print(
    plotExpression(sce.amy, exprs_values = "logcounts", features=genes.top40.t[[i]],
                   x="cellType.split", colour_by="cellType.split", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau20[1:12], length(genes.top40.t[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0(i, " top 40 markers, refined: single-nucleus-level p.w. t-tests, cluster-vs-all"))
  )
  dev.off()
}
