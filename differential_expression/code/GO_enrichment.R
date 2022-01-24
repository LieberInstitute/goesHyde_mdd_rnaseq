library(SummarizedExperiment)
library(jaffelab)
library(org.Hs.eg.db)
library(clusterProfiler)
library(purrr)
library(here)
library(sessioninfo)

#### Load DE results ####
data_type <- c("gene", "exon", "jxn","tx")

files <- here("differential_expression", "data", paste0("qSVA_MDD_",data_type,"_DEresults.rda"))
names(files) <- data_type


allOut <- lapply(files, function(x) get(load(x, verbose = TRUE)))

## sep results only
allOut <- transpose(allOut)
allOut <- allOut$sep

map_depth(allOut, 3, nrow)

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


ont <- c("BP", "CC", "MF", "ALL")
names(ont) <- ont

enrichGO = map(ont, ~compareCluster(geneClusters = gene_sets, 
                    univ = u,
                    OrgDb = "org.Hs.eg.db", 
                    fun = "enrichGO",
                    ont = .x))


# clusterProfiler      * 4.2.0    2021-10-26 [1] Bioconductor
KEGG_DATA <- clusterProfiler:::get_data_from_KEGG_db("hsa")
enrichKEGG = compareCluster(geneClusters = gene_sets, 
                            universe = u,
                            organism = "hsa",
                            fun = "enrichKEGG")

k = enrichKEGG(gene_sets, species = 'hsa')

pdf(file = here("differential_expression", "plots","enrichGO_dotplots.pdf"))
walk2(enrichGO, names(enrichGO), ~dotplot(.x, title = .y))
dev.off()
# 
# save(go, file = here("differential_expression", "data", "go_ontALL.rda"))
# write.csv(go, file = here("differential_expression", "data", "go_ontALL.csv"))


session_info()
