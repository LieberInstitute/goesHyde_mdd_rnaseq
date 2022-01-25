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

#### Functions for extracting gene sets ####
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

my_get_entrez <- function(g){
  e <- unique(ss(g, "\\."))
  entrez <- bitr(e, fromType = "ENSEMBL", toType ="ENTREZID", OrgDb="org.Hs.eg.db")
  return(entrez$ENTREZ)
}

#### Get signif genes####
signif_genes <- map_depth(allOut, 3, ~get_signif(.x, colname = "common_gene_id", return_unique = TRUE))
map_depth(signif_genes, 3, length)

#### Combine All Features ####
## transpose and flatten so features are combined
signif_genes_af <- map(transpose(signif_genes),transpose)
names(signif_genes_af$amyg$MDD)

signif_genes_af_flat <- my_flatten(map_depth(signif_genes_af, 2, unlist))
map_int(signif_genes_af_flat, length)
# amyg.MDD amyg.BPD sacc.MDD sacc.BPD 
# 1141      198     2715      331 

## Extract unique and convert to entrez
gene_sets_af <- map(signif_genes_af_flat, my_get_entrez)

map_int(gene_sets_af, length)
# amyg.MDD amyg.BPD sacc.MDD sacc.BPD 
# 936      138     1668      264 

#### Feature specific for MDD only ####
signif_genes_sf <- transpose(map(signif_genes, transpose))
signif_genes_sf <- signif_genes_sf$MDD
names(signif_genes_sf)

signif_genes_sf_flat <- my_flatten(map_depth(signif_genes_sf, 2, unlist))
map_int(signif_genes_sf_flat, length)
# gene.amyg gene.sacc exon.amyg exon.sacc  jxn.amyg  jxn.sacc   tx.amyg   tx.sacc 
#        50       779       147       973        76       553       868       410 

## Extract unique and convert to entrez
gene_sets_sf <- map(signif_genes_sf_flat, my_get_entrez)

map_int(gene_sets_sf, length)
# gene.amyg gene.sacc exon.amyg exon.sacc  jxn.amyg  jxn.sacc   tx.amyg   tx.sacc 
#        45       708       143       931        76       544       776       364

t(rbind(map_dfr(signif_genes_sf_flat, length), map_dfr(gene_sets_sf, length)))

#### Run Enrichment ####
ont <- c("BP", "CC", "MF", "ALL")
names(ont) <- ont

enrichGO_af = map(ont, ~compareCluster(geneClusters = gene_sets_af, 
                    univ = u,
                    OrgDb = "org.Hs.eg.db", 
                    fun = "enrichGO",
                    ont = .x))

enrichGO_sf = map(ont, ~compareCluster(geneClusters = gene_sets_sf, 
                                       univ = u,
                                       OrgDb = "org.Hs.eg.db", 
                                       fun = "enrichGO",
                                       ont = .x))


# clusterProfiler      * 4.2.0    2021-10-26 [1] Bioconductor
KEGG_DATA <- clusterProfiler:::get_data_from_KEGG_db("hsa")
enrichKEGG = compareCluster(geneClusters = gene_sets_af, 
                            universe = u,
                            organism = "hsa",
                            fun = "enrichKEGG")

k = enrichKEGG(gene_sets_af, species = 'hsa')


#### save Dotplots ####
pdf(file = here("differential_expression", "plots","enrichGO_dotplots.pdf"))
walk2(enrichGO, names(enrichGO), ~dotplot(.x, title = .y))
dev.off()

pdf(file = here("differential_expression", "plots","enrichGO_dotplots_specific.pdf"))
walk2(enrichGO_sf, names(enrichGO_sf), ~dotplot(.x, title = .y))
dev.off()

pdf(file = here("differential_expression", "plots","enrichKEGG_dotplots_all.pdf"))
dotplot(enrichKEGG, title = "enrichKEGG")
dev.off()


walk2(enrichGO, names(enrichGO), function(go_data, name){
  dp <- dotplot(go_data, name)
  ggsave(dp, filename = here("differential_expression", "plots", paste0("enrichGO_dotplots_",name,".png")))
})


# 
save(enrichGO_af, enrichGO_sf, file = here("differential_expression", "data", "enrichGO.rda"))
## Added by Fernando
# write.csv(go, file = here("differential_expression", "data", "go_ontALL.csv"))



session_info()
