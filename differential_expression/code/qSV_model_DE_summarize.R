library(jaffelab)
library(SummarizedExperiment)
library(VennDiagram)
library(tidyverse)
library(sessioninfo)
library(here)
library(EnhancedVolcano)
library(UpSetR)
library(VariantAnnotation)

#### Load Data ####
data_type <- c("gene","exon","jxn","tx")
names(data_type) <- data_type

paths <- map(data_type, ~here("differential_expression","data",paste0("qSVA_MDD_", .x, "_DEresults.rda")))

load(paths$gene, verbose = TRUE)
load(paths$exon, verbose = TRUE)
load(paths$jxn, verbose = TRUE)
load(paths$tx, verbose = TRUE)

allOut <- list(gene = outGene, exon = outExon, jxn = outJxn, tx = outTx)

## Annotate Risk allels ##
load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)
risk_vcf <- readVcf(here("eqtl", "data", "risk_snps", "LIBD_maf01_gwas_MDD.vcf.gz"))
mdd_snps <- read.delim(here("eqtl", "data", "risk_snps", "PGC_depression_genome-wide_significant_makers.txt"))

# risk_gr <- rowRanges(risk_vcf)
risk_gr <- GRanges(seqnames = paste0("chr", mdd_snps$chr), IRanges(mdd_snps$bp, width = 1), feature_id = mdd_snps$markername)
## genes 1KB from risk snp
oo <- findOverlaps(risk_gr, rowRanges(rse_gene), maxgap = 100000)
risk_genes <- unique(rowRanges(rse_gene)[subjectHits(oo),])
length(risk_genes$gencodeID)
##127

#### Summarize Counts ####
get_signif <- function(outFeature, colname = "common_feature_id", cutoff = 0.05, return_unique = FALSE){
      signif <- outFeature[[colname]][outFeature$adj.P.Val < cutoff]
      if(return_unique) signif <- unique(signif)
      signif <- signif[!is.na(signif)]
      return(signif)
}

## Extract lists
signifFeat <- map_depth(allOut, 4, get_signif)
signifGene <- map_depth(allOut, 4, get_signif, colname = "common_gene_id", return_unique = TRUE)
signifFC <- map_depth(allOut, 4, get_signif, colname = "logFC")
signifSymb <- map_depth(allOut, 4, get_signif, colname = "common_gene_symbol")

## Summarize counts
signifFeat_n <- map_depth(signifFeat, 4, length) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(n_features = V1)

signifGene_n <- map_depth(signifGene, 4, length) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(n_genes = V1)

signifRiskGene_n <- map_depth(signifGene, 4, ~sum(.x %in% risk_genes$gencodeID)) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(n_risk_genes = V1)

signifGene_n_up <- map_depth(signifFC, 4, ~sum(.x > 0)) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(Up = V1)

signifGene_n_down <- map_depth(signifFC, 4, ~sum(.x < 0)) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(Down = V1)

signifSymbol <- map_depth(signifSymb, 4, ~paste(unique(.x)[order(unique(.x))], collapse = ", ")) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(Genes = V1)

signifRiskSymbol <- map_depth(signifSymb, 4, ~paste(unique(.x[.x %in% risk_genes$Symbol]), collapse = ", ")) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(RiskGenes = V1)

model_counts <- do.call("cbind",list(signifFeat_n, signifGene_n, signifGene_n_up, signifGene_n_down,signifRiskGene_n, signifRiskSymbol, signifSymbol)) %>%
  rownames_to_column("data") %>%
  separate(data, into = c("Feature", "Model", "Region", "Dx"), sep = "\\.")

model_counts[,1:10]

write_csv(model_counts, here("differential_expression","data","model_counts.csv"))


signifSymb

##### Overlap Between Models ####
my_flatten <- function (x, use.names = TRUE, classes = "ANY") 
{
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

signifFeat_flat <- map(signifFeat,my_flatten)

#### Using map causes bug?
pdf(here("differential_expression","plots",paste0("upset_","gene",".pdf")))
upset(fromList(signifFeat_flat$gene), order.by = "freq", nset = 8)
dev.off()

pdf(here("differential_expression","plots",paste0("upset_exon.pdf")))
upset(fromList(signifFeat_flat$exon), order.by = "freq", nset = 8)
dev.off()

pdf(here("differential_expression","plots",paste0("upset_jxn.pdf")))
upset(fromList(signifFeat_flat$jxn), order.by = "freq", nset = 8)
dev.off()

pdf(here("differential_expression","plots",paste0("upset_tx.pdf")))
upset(fromList(signifFeat_flat$tx), order.by = "freq", nset = 8)
dev.off()

signifGene_flat <- my_flatten(signifGene)
pdf(here("differential_expression","plots",paste0("upset_ALL.pdf")))
upset(fromList(signifGene_flat), order.by = "freq", nset = sum(model_counts$n_genes != 0))
dev.off()

#### Model vs. Model ####
## Compare Model Intersects
signifFeat_model <- map(signifFeat, function(sg) map(transpose(sg), transpose))

get_venn_values <- function(my_sets){
  my_values <- c(length(setdiff(my_sets[[1]], my_sets[[2]])),
                    length(intersect(my_sets[[1]], my_sets[[2]])),
                    length(setdiff(my_sets[[2]], my_sets[[1]])))
  names(my_values) <- c(names(my_sets)[[1]],
                      "Intersect", names(my_sets)[[2]])
  return(my_values)
}

model_venn_values <- as.data.frame(map(signifFeat_model, function(sg) map(sg, ~map(.x, get_venn_values)))) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("cat") %>%
  pivot_longer(!cat)

model_venn_values$name <- factor(model_venn_values$name,  levels = c("sep", "Intersect", "sep_cf"))

model_venn_tiles <- ggplot(model_venn_values, aes(x = name, y = cat, fill = value)) +
  geom_tile(color = "lightgrey") +
  geom_text(aes(label = value), color = "white") +
  labs(x = "n features", y = NULL)
  
ggsave(model_venn_tiles, filename = here("differential_expression","plots", "model_venn_tiles.png"))

## Graph t-stats
allOut_model <- map(allOut, function(x) map(transpose(x), transpose))
names(allOut_model$gene$amyg$MDD)

head(do.call("cbind", allOut_model$gene$amyg$MDD))

allOut_model_compare <- map_depth(allOut_model, 3, do.call, what = "cbind")

sig_colors <- c(RColorBrewer::brewer.pal(3, "Set1"),"dark grey")
names(sig_colors) <- c("sig_Both", "sig_no-deconvo", "sig_deconvo", "None")

model_compare_plot <- function(model_df, title){

  t_cor <- model_df %>%
    summarize(cor = cor(sep.t, sep_cf.t, method = "spearman")) %>%
    mutate(cor_anno = paste0("rho==", format(round(cor, 2), nsmall = 2)))
  
  model_df <- model_df %>% 
    mutate(Signif = case_when(sep.adj.P.Val < 0.05 & sep_cf.adj.P.Val < 0.05 ~"sig_Both",
                              sep.adj.P.Val < 0.05 ~ "sig_no-deconvo",
                              sep_cf.adj.P.Val < 0.05 ~ "sig_deconvo",
                              TRUE ~ "None"))
  
  model_plot <- ggplot(model_df, aes(x = sep.t, y = sep_cf.t, color = Signif)) +
    geom_point(alpha = 0.5, size = 0.5) +
    labs(x = "t-stats DE Model", y = "t-stats DE Model with Cell Fractions", 
         title = title, subtitle = t_cor$cor_anno) +
    scale_color_manual(values = sig_colors)+
    theme_bw()
  
  message(title)  
  print(model_plot)
}

pdf(here("differential_expression","plots","qSV_model_scatter.pdf"))
walk2(allOut_model_compare, names(allOut_model_compare), function(feat_df, feat_name){
  walk2(feat_df, names(feat_df), function(region_df, region_name){
    walk2(region_df, names(region_df), ~model_compare_plot(.x, title = paste(feat_name, region_name,.y,"vs. Control")))
  })
})
dev.off()

## Compare Dx or Region?

#### Volcano Plots ####

pdf(here("differential_expression","plots","qSV_model_volcano.pdf"))
map2(allOut, c("Gene", "Exon", "Jxn", "Tx"), function(featOut, featName){
  map2(featOut, names(featOut), function(modelOut, modelName){
    map2(modelOut, names(modelOut), function(regionOut, regionName){
      map2(regionOut, names(regionOut),
           ~EnhancedVolcano(.x,
                            pCutoff = 0.05,
                            y = "adj.P.Val",
                            x = "logFC",
                            lab = .x$Symbol,
                            title = paste(featName,"-", regionName),
                            subtitle = paste(modelName, "model,", .y,"vs. Control"))
      )
    })
  })
})

dev.off()

#### Boxplots ####

# ranks_to_check <- c(1:50, 201:250, 401:450, 601:650, 801:850, 1001:1050)
ranks_to_check <- c(1:50)
length(ranks_to_check)

outGene_test <- outGene$sep$amyg$MDD %>%
  arrange(adj.P.Val) %>%
  mutate(row_num = row_number())

