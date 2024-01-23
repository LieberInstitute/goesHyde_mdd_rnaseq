library(SummarizedExperiment)
library(jaffelab)
library(org.Hs.eg.db)
library(clusterProfiler)
library(purrr)
library(here)
library(sessioninfo)

#here() starts at /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq

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


          [,1] [,2]
gene.amyg   50   45
gene.sacc  779  708
exon.amyg  147  143
exon.sacc  973  931
jxn.amyg    76   76
jxn.sacc   553  544
tx.amyg    868  776
tx.sacc    410  364


 str(gene_sets_af)
List of 4
 $ amyg.MDD: chr [1:936] "7398" "50999" "65005" "6283" ...
 $ amyg.BPD: chr [1:138] "5351" "163732" "114625" "6513" ...
 $ sacc.MDD: chr [1:1668] "653635" "57801" "100288175" "142678" ...
 $ sacc.BPD: chr [1:264] "1325" "55920" "126917" "440574" ...
>

#### Run Enrichment ####
ont <- c("BP", "CC", "MF", "ALL")
names(ont) <- ont

enrichGO_af = map(ont, ~compareCluster(geneClusters = gene_sets_af, 
                    univ = u,
                    OrgDb = "org.Hs.eg.db", 
                    fun = "enrichGO",
                    ont = .x))


#Per feature
enrichGO_sf = map(ont, ~compareCluster(geneClusters = gene_sets_sf, 
                                       univ = u,
                                       OrgDb = "org.Hs.eg.db", 
                                       fun = "enrichGO",
                                       ont = .x))

#enrichGO_sf has eight levels -- each feature for MDD and BD

#use enrichGO_af for combined across features

enrichGO_af_simplify = map(enrichGO_af, ~clusterProfiler::simplify((.x), cutoff=0.5, by="qvalue", select_fun=min))

enrichGO_sf_simplify = map(enrichGO_sf, ~clusterProfiler::simplify((.x), cutoff=0.5, by="qvalue", select_fun=min))

write.csv(enrichGO_sf@compareClusterResult, file = here("differential_expression", "data", "go_ontALL_020521.csv"))


emapplot.compareClusterResult(enrichGO_af_simplify)

edo = pairwise_termsim(enrichGO_af_simplify$BP)

emapplot(edo, cex_category=0.6, cex_label_category = 0.5, showCategory = 40)





#REVIGO RESULTS -----

sacc_revigo4 <- c("GO:0002181", "GO:0060627", "GO:0070997", "GO:0016358", "GO:0099504", "GO:0030010", "GO:0043112", "GO:0061684", "GO:0043393", "GO:1902414", "GO:0001505", "GO:0090257", "GO:0009205", "GO:0006836", "GO:0043484", "GO:0060401", "GO:0050804", "GO:0051235", "GO:1901617", "GO:0031503", "GO:0043254", "GO:0007264", "GO:0032970", "GO:0022604", "GO:0018107", "GO:0034250", "GO:0035304", "GO:0048024", "GO:0010038", "GO:0010644", "GO:1903311", "GO:0051258", "GO:0008380", "GO:0006470", "GO:0050808", "GO:0034248", "GO:0006684", "GO:0016311", "GO:0097242", "GO:1902570", "GO:0022411", "GO:0005844", "GO:0022626", "GO:0043209", "GO:0150034", "GO:0098984", "GO:0030427", "GO:0043025", "GO:0031252", "GO:0098562", "GO:0042383", "GO:0030139", "GO:0097038", "GO:0005938", "GO:0045259", "GO:0036464", "GO:0070603", "GO:0003735", "GO:0003924", "GO:0015662", "GO:0030695", "GO:0048156", "GO:0003746", "GO:0003697", "GO:0051020", "GO:0070325", "GO:0044325", "GO:0045296", "GO:0017124")
sacc_revigo4_BP <- c("GO:0002181", "GO:0060627", "GO:0070997", "GO:0016358", "GO:0099504", "GO:0030010", "GO:0043112", "GO:0061684", "GO:0043393", "GO:1902414", "GO:0001505", "GO:0090257", "GO:0009205", "GO:0006836", "GO:0043484", "GO:0060401", "GO:0050804", "GO:0051235", "GO:1901617", "GO:0031503", "GO:0043254", "GO:0007264", "GO:0032970", "GO:0022604", "GO:0018107", "GO:0034250", "GO:0035304", "GO:0048024", "GO:0010038", "GO:0010644", "GO:1903311", "GO:0051258", "GO:0008380", "GO:0006470", "GO:0050808", "GO:0034248", "GO:0006684", "GO:0016311", "GO:0097242", "GO:1902570", "GO:0022411") 
sacc_revigo4_CC <- c("GO:0005844", "GO:0022626", "GO:0043209", "GO:0150034", "GO:0098984", "GO:0030427", "GO:0043025", "GO:0031252", "GO:0098562", "GO:0042383", "GO:0030139", "GO:0097038", "GO:0005938", "GO:0045259", "GO:0036464", "GO:0070603")
sacc_revigo4_MF <- c("GO:0003735", "GO:0003924", "GO:0015662", "GO:0030695", "GO:0048156", "GO:0003746", "GO:0003697", "GO:0051020", "GO:0070325", "GO:0044325", "GO:0045296", "GO:0017124")


amyg_revigo1_BP <- c("GO:0043484", "GO:0070936", "GO:0051258", "GO:0008380", "GO:0043161", "GO:0000381", "GO:0042176", "GO:1903311", "GO:0071044", "GO:0000209", "GO:0000398", "GO:0051865", "GO:0000380", "GO:0010498", "GO:0006397") 
amyg_revigo1_CC <- c("GO:0150034", "GO:0030427", "GO:0014069", "GO:0022626", "GO:0005681", "GO:0035371", "GO:0071013", "GO:0098984", "GO:0044309", "GO:0043197", "GO:0099572", "GO:0043025", "GO:0030426")
amyg_revigo1_MF <- c("GO:0051010", "GO:0061631", "GO:0019787", "GO:0031625", "GO:0061650", "GO:0044389", "GO:0004842")

amyg_revigo2_BP <- c("GO:0043484", "GO:0070936", "GO:0051258", "GO:0008380", "GO:0043161", "GO:0000381", "GO:0042176", "GO:1903311", "GO:0071044", "GO:0000209", "GO:0000398", "GO:0051865")
amyg_revigo2_CC <- c("GO:0150034", "GO:0030427", "GO:0014069", "GO:0022626", "GO:0005681", "GO:0035371", "GO:0071013", "GO:0098984")
amyg_revigo2_MF <- c("GO:0051010", "GO:0061631", "GO:0019787", "GO:0031625", "GO:0061650")
#--------

#USE GG PLOT INTEAD OF CLUSTER PROFILER-----
d = as.data.frame(enrichGO_af$ALL@compareClusterResult)

table(d$Cluster)

## make the ratio numeric
d$GeneProp = as.numeric(ss(d$GeneRatio, "/"))/as.numeric(ss(d$GeneRatio, "/",2))

str(d)
## you can filter way easier
d[d=="proteasome−mediated ubiquitin−dependent protein catabolic process"] <- "proteasome/ubiquitin process"
d[d=="proton−transporting ATP synthase complex"] <- "proton ATP synthase complex"

saccBP <- d %>% filter(Cluster == "sacc.MDD") %>% filter(ID %in% sacc_revigo4_BP) 
saccCC <- d %>% filter(Cluster == "sacc.MDD") %>% filter(ID %in% sacc_revigo4_CC) 
saccMF <- d %>% filter(Cluster == "sacc.MDD") %>% filter(ID %in% sacc_revigo4_MF) 

amygBP <- d %>% filter(Cluster == "amyg.MDD") %>% filter(ID %in% amyg_revigo1_BP) 
amygCC <- d %>% filter(Cluster == "amyg.MDD") %>% filter(ID %in% amyg_revigo1_CC) 
amygMF <- d %>% filter(Cluster == "amyg.MDD") %>% filter(ID %in% amyg_revigo1_MF) 

sacc <- d %>% filter(Cluster == "sacc.MDD") %>% filter(ID %in% sacc_revigo4_BP) 

temp = sacc_revigo4_BP + sacc_revigo4_CC + sacc_revigo4_MF

## and plot


g_saccBP <- ggplot(saccBP, aes(x=-log10(qvalue), y=reorder(Description, -qvalue))) +  geom_point(aes(size = GeneProp), colour = "red", alpha = 0.9, shape = 20) + theme_dose() + ylab(NULL) + xlab(label = '-log FDR')  + theme(axis.text.y=element_text(size =8)) + ggtitle("MDD enrichment GO Biological Processes in sACC")
g_saccCC <- ggplot(saccCC, aes(x=-log10(qvalue), y=reorder(Description, -qvalue))) +  geom_point(aes(size = GeneProp), colour = "red", alpha = 0.9, shape = 20) + theme_dose() + ylab(NULL) + xlab(label = '-log FDR')  + ggtitle("MDD enrichment GO Cellular Component in sACC")
g_saccMF <- ggplot(saccMF, aes(x=-log10(qvalue), y=reorder(Description, -qvalue))) +  geom_point(aes(size = GeneProp), colour = "red", alpha = 0.9, shape = 20) + theme_dose() + ylab(NULL) + xlab(label = '-log FDR')  + ggtitle("MDD enrichment GO Molecular Function in sACC")

g_amygBP <- ggplot(amygBP, aes(x=-log10(qvalue), y=reorder(Description, -qvalue))) +  geom_point(aes(size = GeneProp), colour = "red", alpha = 0.9, shape = 20) + theme_dose() + ylab(NULL) + xlab(label = '-log FDR')   + ggtitle("MDD enrichment GO Biological Processes in Amygdala")
g_amygCC <- ggplot(amygCC, aes(x=-log10(qvalue), y=reorder(Description, -qvalue))) +  geom_point(aes(size = GeneProp), colour = "red", alpha = 0.9, shape = 20) + theme_dose() + ylab(NULL) + xlab(label = '-log FDR')   + ggtitle("MDD enrichment GO Cellular Component in Amygdala")
g_amygMF <- ggplot(amygMF, aes(x=-log10(qvalue), y=reorder(Description, -qvalue))) +  geom_point(aes(size = GeneProp), colour = "red", alpha = 0.9, shape = 20) + theme_dose() + ylab(NULL) + xlab(label = '-log FDR')   + ggtitle("MDD enrichment GO Molecular Function in Amygdala")

coord_fixed()

#------

library("tidyverse")
library(grid)
library(cowplot)

#     xlab("Modules") + ylab("t-statistic") + ggtitle("MDD Associated Coexpression Modules in the sACC") +

pdf(file = here("differential_expression", "plots","combined_features_GOenrich_020622.pdf"),  width = 10, height = 6)
g_saccBP 
g_saccCC  + xlim(0, 8)
g_saccMF  + xlim(0, 5)

g_amygBP + xlim(0, 3)
g_amygCC + xlim(0, 3)
g_amygMF + xlim(0, 3)
dev.off()
#same.size.ggplot(c("g_saccBP", "g_saccCC", "g_saccMF"))

#grid.draw(rbind(ggplotGrob(g_saccBP), ggplotGrob(g_saccCC), ggplotGrob(g_saccMF)))
pdf(file = here("differential_expression", "plots","combined_features_GOenrich_020622.pdf"),  width = 10, height = 6)

plot_grid(g_saccBP, g_saccCC, g_saccMF, align = "h", nrow = 3, rel_heights = c(1,1,1))
plot_grid(g_saccCC, align = "h", nrow = 1, rel_heights = c(1/4, 1/4, 1/4))

dev.off()

pdf(file = here("differential_expression", "plots","test.pdf"),  width = 10, height = 6)
#ggplot(saccBP, aes(x=-log10(qvalue), y=reorder(Description, -qvalue))) +  geom_bar(stat ="identity", color = saccBP$qvalue) + scale_colour_gradient(limits=c(0, max(d$qvalue)), low="red",high="blue") + theme_dose() + ylab(NULL) + theme(axis.text.y=element_text(size =8)) + ggtitle("MDD enrichment GO Biological Processes in sACC")


ggplot(saccMF, aes(x=-log10(qvalue), y=reorder(Description, -qvalue))) + 
    geom_point(aes(size = GeneProp), colour = "red", alpha = 0.5, shape = 20) + 
#    scale_colour_gradient(limits=c(0, max(d$qvalue)), low="red",high="blue") +
    theme_dose() + ylab(NULL) + xlab(label = '-log FDR') +
    theme(axis.text.y=element_text(size =8)) + 
    ggtitle("MDD enrichment GO Biological Processes in sACC") +
    xlim(0, 3)

dev.off()


#### Run Enrichment ####




# clusterProfiler      * 4.2.0    2021-10-26 [1] Bioconductor
KEGG_DATA <- clusterProfiler:::get_data_from_KEGG_db("hsa")
enrichKEGG = compareCluster(geneClusters = gene_sets_af, 
                            universe = u,
                            organism = "hsa",
                            fun = "enrichKEGG")

k = enrichKEGG(gene_sets_af, species = 'hsa')


 
save(enrichGO_af, enrichGO_sf, file = here("differential_expression", "data", "enrichGO.rda"))
## Added by Fernando
# write.csv(go, file = here("differential_expression", "data", "go_ontALL.csv"))

GOenrich <- load(here("differential_expression", "data", "enrichGO.rda"))
str(enrichGO_sf)
names(enrichGO_sf)
names(enrichGO_af)


session_info()
