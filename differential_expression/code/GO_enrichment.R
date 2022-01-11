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

#load objects
load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)
load(here('exprs_cutoff','rse_exon.Rdata'), verbose=TRUE)
load(here('exprs_cutoff','rse_jxn.Rdata'), verbose=TRUE)
load(here('exprs_cutoff','rse_tx.Rdata'), verbose=TRUE)

#### Build Expression Map ####
## modify rowRanges
rowRanges(rse_gene)$Type <- "Gene"
rowRanges(rse_exon)$Type <- "Exon"

# jxn "newGeneID"     "newGeneSymbol"
colnames(mcols(rowRanges(rse_jxn)))[13:14] <- c("gencodeID", "Symbol")
rowRanges(rse_jxn)$Type <- "Junction"

# Exon "gene_id"   "gene_name"
colnames(mcols(rowRanges(rse_tx)))[c(5,8)] <- c("gencodeID", "Symbol")
rowRanges(rse_tx)$Class <- "InGen"
rowRanges(rse_tx)$Type <- "Transcript"

## Expression Map
name <- c("gencodeID", "Symbol", "Type", "Class")
exprsMap <- rbind(as.data.frame(rowRanges(rse_gene))[,name],
                 as.data.frame(rowRanges(rse_exon))[,name],
                 as.data.frame(rowRanges(rse_jxn))[,name],
                 as.data.frame(rowRanges(rse_tx))[,name])
exprsMap <- DataFrame(exprsMap)

nrow(exprsMap)
# [1] 810029

#### Load DE results ####
data_type <- c("gene", "exon", "jxn","tx")

files <- here("differential_expression", "data", paste0("qSVA_MDD_",data_type,"_DEresults.rda"))
names(files) <- data_type


allOut <- lapply(files, function(x) get(load(x, verbose = TRUE)))
map_depth(allOut, 4, ~"EntrezID" %in% colnames(.x))

allEntrez <- map_depth(allOut, 4, ~.x$EntrezID)
allEntrez_ul <- unlist(allEntrez)
u <- as.character(allEntrez_ul[!is.na(allEntrez_ul)])
length(u)

sigOut <- map_depth(allOut, 4, ~.x[.x$adj.P.Val < 0.05,])
map_depth(sigOut, 4, nrow)

sigGenes <- map_depth(sigOut, 4, ~.x$EntrezID)

g_test <- map(sigGenes, ~pluck(.x, "sep","sacc","MDD"))
g_test2 <- map(g_test, ~as.character(unique(.x[!is.na(.x)])))
lengths(g_test2[1:2])
# gene exon 
# 678  906
is.character(g_test2[[1]])
is.character(g_test2)

go = compareCluster(geneClusters = g_test2[[1]], 
                    univ = u,
                    OrgDb = "org.Hs.eg.db", 
                    ont = "ALL",
                    readable = TRUE, 
                    pvalueCutoff = 1,
                    qvalueCutoff = 1)

# Error in (function (cl, name, valueClass)  : 
#            assignment of an object of class “NULL” is not valid for @‘fun’ in an object of class “compareClusterResult”; is(value, "character") is not TRUE
