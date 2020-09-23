## Load required R packages
library("BiocFileCache")
library("S4Vectors")
library("ggplot2")
library("plotly")
library("sessioninfo")
library("SummarizedExperiment")
library("dplyr")
library("reshape2")
library("rlang")

plotxy <- function(xvar, yvar, df = pd_m, color_var = "BrainRegion", shape_var = "Experiment") {
  ggplot(df, aes(x = !!sym(xvar), y = !!sym(yvar), 
                 color = !!sym(color_var), shape = !!sym(shape_var))) + 
    geom_point() + 
    labs(x = xvar, y = yvar) 
}

#load data 
qcresults <- read.csv("qc_dropping_results.csv",row.names=1)

load("../data/rse_gene_raw_GoesZandi.rda", verbose = TRUE)
pd <- colData(rse_gene)

id_cols <- c("RNum","BrNum","Experiment","BrainRegion")
metric_cols <- c("RIN", "overallMapRate", "totalAssignedGene", "mitoRate", "numReads" ,"rRNA_rate","Plate", "ERCCsumLogErr")

pd_m <- pd %>% 
  as.data.frame() %>%
  select(-BrainRegion) %>%
  left_join(qcresults) %>%
  select(all_of(c(id_cols, metric_cols))) %>%
  mutate(log10_numReads = log10(numReads))

## ERCC plotly
y_metrics <- c("overallMapRate", "totalAssignedGene", "log10_numReads","RIN","mitoRate", "rRNA_rate")
x_metrics <- "ERCCsumLogErr"

pd_m_key <- highlight_key(pd_m, ~ RNum)
ercc_vs_plots<- mapply(plotxy, x_metrics, y_metrics, MoreArgs = list(df = pd_m_key), SIMPLIFY = FALSE)

# add existing cutoffs
ercc_vs_plots[[1]] <- ercc_vs_plots[[1]] + geom_hline(yintercept=0.5, linetype="dashed")
ercc_vs_plots[[2]] <- ercc_vs_plots[[2]] + geom_hline(yintercept=0.3, linetype="dashed")
ercc_vs_plots[[3]] <- ercc_vs_plots[[3]] + geom_hline(yintercept=7.25, linetype="dashed")
ercc_vs_plots[[6]] <- ercc_vs_plots[[6]] + geom_hline(yintercept=1e-3, linetype="dashed")

ercc_vs_plots <- lapply(ercc_vs_plots, function(x) x +theme(legend.position= "none"))

p_ercc_vs_plots <- lapply(ercc_vs_plots, ggplotly)
p_merged <- subplot(p_ercc_vs_plots, 
                    margin = 0.05,
                    titleY = TRUE,
                    titleX = TRUE,
                    shareX = TRUE,
                    nrows = 2)

htmlwidgets::saveWidget(highlight(p_merged,
                                  on = "plotly_click",
                                  off = "plotly_doubleclick",
                                  selectize = TRUE,
                                  dynamic = TRUE,
                                  persistent = FALSE),
                        file = 'ERCCsumLogErr_vs_metrics.html')

# sgejobs::job_single('qc_plotly', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript qc_plotly.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
