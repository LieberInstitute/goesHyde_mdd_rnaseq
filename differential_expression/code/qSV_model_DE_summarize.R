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

#### Summarize Counts ####
#model, region, Dx
summarize_DE <- function(outFeature, feature_name){
  feature_df <- map_dfr(outFeature, function(outModel){
    model_df <- map_dfr(outModel, function(outRegion){
      region_df <- map_dfr(outRegion, ~sum(.x$adj.P.Val < 0.05, na.rm = TRUE))
    })
    model_df$region <- names(outModel)
    return(model_df)
  })
  feature_df$model <- rep(names(outFeature), each = 2)
  feature_df <- feature_df %>%
    pivot_longer(!c(region, model), names_to = "Dx", values_to =  feature_name)
  
  return(feature_df)
  }

model_counts <- summarize_DE(outGene, "Gene")  %>%
  left_join(summarize_DE(outExon, "Exon"))%>%
  left_join(summarize_DE(outJxn, "Jxn"))%>%
  left_join(summarize_DE(outTx, "Tx"))

# region model  Dx     Gene  Exon   Jxn    Tx
# <chr>  <chr>  <chr> <int> <int> <int> <int>
# 1 amyg   sep    MDD       3     4    47    66
# 2 amyg   sep    BPD      82   159    23     7
# 3 sacc   sep    MDD     392   430   470    43
# 4 sacc   sep    BPD     326   119    34    18
# 5 amyg   sep_cf MDD       4     0    47     9
# 6 amyg   sep_cf BPD      12    44    16     0
# 7 sacc   sep_cf MDD      95    20   409    26
# 8 sacc   sep_cf BPD      21     5    37     1
 
write_csv(model_counts, here("differential_expression","data","model_counts.csv"))

##### Overlap Between Models ####
map(list(outGene, outExon, outJxn, outTx), ~head(.x$sep$amyg$MDD))

get_signif_genes <- function(outFeature){
  feature_genes <- map(outFeature, function(outModel){
    model_genes <- map(outModel, function(outRegion){
      region_genes <- map(outRegion, ~.x$common_gene_id[.x$adj.P.Val < 0.05])
    })
  })
  return(feature_genes)
}
signifGene <- get_signif_genes(outGene)

flatten(flatten(signifGene))

map(signifGene$sep, ~length(intersect(.x$MDD, .x$BPD)))
map(signifGene,function(mod) map(mod, ~length(intersect(.x$MDD, .x$BPD))))

#### Volcano Plots ####

pdf(here("differential_expression","plots","qSV_model_volcano.pdf"))
map2(list(outGene, outExon, outJxn, outTx), c("Gene", "Exon", "Jxn", "Tx"), function(featOut, featName){
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






