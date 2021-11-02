library(jaffelab)
library(SummarizedExperiment)
library(VennDiagram)
library(tidyverse)
library(sessioninfo)
library(here)
library(EnhancedVolcano)
library(UpSetR)

#### Load Data ####
data_type <- c("gene","exon","jxn","tx")
names(data_type) <- data_type

paths <- map(data_type, ~here("differential_expression","data",paste0("qSVA_MDD_", .x, "_DEresults.rda")))

load(paths$gene, verbose = TRUE)
load(paths$exon, verbose = TRUE)
load(paths$jxn, verbose = TRUE)
load(paths$tx, verbose = TRUE)

allOut <- list(gene = outGene, exon = outExon, jxn = outJxn, tx = outTx)

#### Summarize Counts ####
get_signif <- function(outFeature, colname = "common_feature_id", cutoff = 0.05){
      signif <- unique(outFeature[[colname]][outFeature$adj.P.Val < cutoff])
      signif <- signif[!is.na(signif)]
      return(signif)
}

## Extract lists
signifFeat <- map_depth(allOut, 4, get_signif)
signifGene <- map_depth(allOut, 4, get_signif, colname = "common_gene_id")
signifFC <- map_depth(allOut, 4, get_signif, colname = "logFC")

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

model_counts <- do.call("cbind",list(signifFeat_n, signifGene_n, signifGene_n_up, signifGene_n_down)) %>%
  rownames_to_column("data") %>%
  separate(data, into = c("Feature", "Model", "Region", "Dx"), sep = "\\.")

write_csv(model_counts, here("differential_expression","data","model_counts.csv"))

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

# pdf(here("differential_expression","plots",paste0("upset_tx.pdf")))
# upset(fromList(signifFeat_flat$tx), order.by = "freq")
# dev.off()

signifGene_flat <- my_flatten(signifGene)
pdf(here("differential_expression","plots",paste0("upset_ALL.pdf")))
upset(fromList(signifGene_flat), order.by = "freq", nset = 24)
dev.off()

## Compare Model Intersects
signifGene_model <- map(signifGene, function(sg) map(transpose(sg), transpose))

get_venn_values <- function(my_sets){
  my_values <- c(length(setdiff(my_sets[[1]], my_sets[[2]])),
                    length(intersect(my_sets[[1]], my_sets[[2]])),
                    length(setdiff(my_sets[[2]], my_sets[[1]])))
  names(my_values) <- c(names(my_sets)[[1]],
                      "Intersect", names(my_sets)[[2]])
  return(my_values)
}

model_venn_values <- as.data.frame(map(signifGene_model, function(sg) map(sg, ~map(.x, get_venn_values)))) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("cat") %>%
  pivot_longer(!cat)

model_venn_values$name <- factor(model_venn_values$name,  levels = c("sep", "Intersect", "sep_cf"))

model_venn_tiles <- ggplot(model_venn_values, aes(x = name, y = cat, fill = value)) +
  geom_tile(color = "lightgrey") +
  geom_text(aes(label = value), color = "white")
  
ggsave(model_venn_tiles, filename = here("differential_expression","plots", "model_venn_tiles.png"))

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






