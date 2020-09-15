## Load required R packages
library("BiocFileCache")
library("S4Vectors")
library("ggplot2")
library("plotly")
library("sessioninfo")
library("SummarizedExperiment")
library("dplyr")
library("reshape2")

#load data 
qcresults <- read.csv("qc_dropping_results.csv",row.names=1)

load("../data/rse_gene_raw_GoesZandi.rda", verbose = TRUE)
pd <- colData(rse_gene)

id_cols <- c("RNum","BrNum","Experiment","BrainRegion")
metric_cols <- c("RIN", "overallMapRate", "totalAssignedGene", "mitoRate", "rRNA_rate","Plate", "ERCCsumLogErr")

pd_m <- pd %>% 
  as.data.frame() %>%
  left_join(qcresults) %>%
  select(all_of(c(id_cols, metric_cols)))


## Let's make a "highlighted" table
pd_m$key <- pd_m$RNum
pd_key <- highlight_key(pd_m, ~ key)

## Plot mitoRate vs. AssignedGene
gg_mito_vs_gene <-
  ggplot(pd_key,
         aes(x = mitoRate, y = totalAssignedGene, color = BrainRegion)) + geom_point()

## Make a second plot
gg_mito_vs_RIN <-
  ggplot(pd_key, aes(x = mitoRate, y = RIN, color = BrainRegion)) + geom_point()

## Convert them to interactive plots
p_mito_vs_gene <- ggplotly(gg_mito_vs_gene)
p_mito_vs_RIN <- ggplotly(gg_mito_vs_RIN)

## Now group them together
p_merged <- subplot(
  p_mean_mito_vs_mean_gene,
  p_mean_mito_vs_mean_RIN,
  nrows = 1,
  shareX = TRUE,
  shareY = FALSE,
  which_layout = 2
)

## Save the version you liked
htmlwidgets::saveWidget(highlight(p_merged, on = "plotly_click", off = "plotly_doubleclick"),
                        file = 'index.html')

