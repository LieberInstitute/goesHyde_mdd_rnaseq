library(SummarizedExperiment)
library(edgeR)
library(limma)
library(recount)
library(jaffelab)
library(IRanges)
library(org.Hs.eg.db)
library(clusterProfiler)
library(purrr)
library(here)

# #load objects
# load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)
# load(here('exprs_cutoff','rse_exon.Rdata'), verbose=TRUE)
# load(here('exprs_cutoff','rse_jxn.Rdata'), verbose=TRUE)
# load(here('exprs_cutoff','rse_tx.Rdata'), verbose=TRUE)
# 
# #### Build Expression Map ####
# ## modify rowRanges
# rowRanges(rse_gene)$Type <- "Gene"
# rowRanges(rse_exon)$Type <- "Exon"
# 
# # jxn "newGeneID"     "newGeneSymbol"
# colnames(mcols(rowRanges(rse_jxn)))[13:14] <- c("gencodeID", "Symbol")
# rowRanges(rse_jxn)$Type <- "Junction"
# 
# # Exon "gene_id"   "gene_name"
# colnames(mcols(rowRanges(rse_tx)))[c(5,8)] <- c("gencodeID", "Symbol")
# rowRanges(rse_tx)$Class <- "InGen"
# rowRanges(rse_tx)$Type <- "Transcript"
# 
# ## Expression Map
# name <- c("gencodeID", "Symbol", "Type", "Class")
# exprsMap <- rbind(as.data.frame(rowRanges(rse_gene))[,name],
#                  as.data.frame(rowRanges(rse_exon))[,name],
#                  as.data.frame(rowRanges(rse_jxn))[,name],
#                  as.data.frame(rowRanges(rse_tx))[,name])
# exprsMap <- DataFrame(exprsMap)
# 
# nrow(exprsMap)
# # [1] 810029


#### Load DE results ####
data_type <- c("gene", "exon", "jxn","tx")

files <- here("differential_expression", "data", paste0("qSVA_MDD_",data_type,"_DEresults.rda"))
names(files) <- data_type


allOut <- lapply(files, function(x) get(load(x, verbose = TRUE)))

## sep results only
allOut <- transpose(allOut)
allOut <- allOut$sep

map_depth(allOut, 3, nrow)
map_depth(allOut, 3, colnames)

## Define Universe
all_gencode <- map_depth(allOut, 3, "common_gene_id")
map_depth(all_gencode, 3, length)
head(all_gencode$gene$amyg$MDD)

length(unlist(all_gencode))
# [1] 3240116

## all ENSEMBL
all_ensembl <- unique(ss(unlist(all_gencode),"\\."))
length(all_ensembl)
# [1] 34635

all_entrez <- bitr(all_ensembl, fromType = "ENSEMBL", toType ="ENTREZID", OrgDb="org.Hs.eg.db")
# 32.71% of input gene IDs are fail to map...
nrow(all_entrez)
# [1] 23468
u <- all_entrez$ENTREZ

## Extract gene sets
get_signif <- function(outFeature, colname = "common_feature_id", cutoff = 0.05, return_unique = FALSE){
  signif <- outFeature[[colname]][outFeature$adj.P.Val < cutoff]
  if(return_unique) signif <- unique(signif)
  signif <- signif[!is.na(signif)]
  return(signif)
}
my_flatten <- function (x, use.names = TRUE, classes = "ANY") {
  #' Source taken from rlist::list.flatten
  len <- sum(rapply(x, function(x) 1L, classes = classes))
  y <- vector("list", len)
  i <- 0L
  items <- rapply(x, function(x) {
    i <<- i + 1L
    y[[i]] <<- x
    TRUE
  }, classes = classes)
  if (use.names && !is.null(nm <- names(items))) 
    names(y) <- nm
  y
}

signif_genes <- map_depth(allOut, 3, ~get_signif(.x, colname = "common_gene_id", return_unique = TRUE))
map_depth(signif_genes, 3, length)

## transpose and flatten so features are combined
signif_genes2 <- map(transpose(signif_genes),transpose)
names(signif_genes2$amyg$MDD)

signif_genes3 <- my_flatten(map_depth(signif_genes2, 2, unlist))
map_int(signif_genes3, length)
# amyg.MDD amyg.BPD sacc.MDD sacc.BPD 
# 1141      198     2715      331 

## Extract unique and convert to entrez
gene_sets <- map(signif_genes3, function(g){
  e <- unique(ss(g, "\\."))
  entrez <- bitr(e, fromType = "ENSEMBL", toType ="ENTREZID", OrgDb="org.Hs.eg.db")
  return(entrez$ENTREZ)
})

map_int(gene_sets, length)
# amyg.MDD amyg.BPD sacc.MDD sacc.BPD 
# 936      138     1668      264 

map(gene_sets, class)

ck <- compareCluster(geneCluster = gene_sets, 
                     fun = enrichGO, 
                     OrgDb = "org.Hs.eg.db", univ = u)
dotplot(ck)

go = compareCluster(geneClusters = gene_sets, 
                    univ = u,
                    OrgDb = "org.Hs.eg.db", 
                    # ont = "ALL",
                    # readable = TRUE, 
                    # pvalueCutoff = 1,
                    # qvalueCutoff = 1,
                    fun = enrichGO)

dotplot(go)

save(go, file = here("differential_expression", "data", "go.rda"))
