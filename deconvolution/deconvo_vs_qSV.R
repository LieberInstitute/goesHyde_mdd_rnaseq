library(SummarizedExperiment)
library(jaffelab)
library(here)
library(reshape2)
library(purrr)
library(tidyverse)
library(broom)
library(viridis)
library(sessioninfo)
library(RColorBrewer)

#### Load Data ####
load(here("deconvolution","data","est_prop_top5_long.Rdata"), verbose = TRUE)
load(here("differential_expression" ,"qSV_mat.Rdata"), verbose = TRUE)

qSV_long <- melt(qSV_mat) %>% rename(sample = Var1, qSV = Var2, qSV_value = value)

## Bind with qSV table
map_int(est_prop_long, nrow)
# sacc_broad    amyg_broad sacc_specific amyg_specific 
# 3306          3240          5510          4860
est_prop_qsv <- map(est_prop_long, ~.x %>% 
                      select(sample, cell_type, prop) %>%
                      left_join(qSV_long, by = "sample"))


#### Calculate p-values ####

prop_qSV_fit <- map(est_prop_qsv,~.x %>% group_by(cell_type, qSV) %>%
                      do(fitQSV = tidy(lm(prop ~ qSV_value-1, data = .))) %>% 
                      unnest(fitQSV) %>%
                      mutate(p.bonf = p.adjust(p.value, "bonf"),
                             p.bonf.sig = p.bonf < 0.05,
                             p.bonf.cat = cut(p.bonf, 
                                              breaks = c(1,0.05, 0.01, 0.005, 0),
                                              labels = c("< 0.005","< 0.01", "< 0.05", "> 0.05")),
                             p.fdr = p.adjust(p.value, "fdr"))
)

map_int(prop_qSV_fit, nrow)
# sacc_broad    amyg_broad sacc_specific amyg_specific 
# 156           156           260           234 
map_int(prop_qSV_fit, ~sum(.x$p.value < 0.05))
# sacc_broad    amyg_broad sacc_specific amyg_specific 
# 67            81           119           114 
map_int(prop_qSV_fit, ~sum(.x$p.bonf.sig))
# sacc_broad    amyg_broad sacc_specific amyg_specific 
# 44            44            68            62 
map_int(prop_qSV_fit, ~sum(.x$p.fdr < 0.05))
# sacc_broad    amyg_broad sacc_specific amyg_specific 
# 60            70           105           102 

#### Save results ###
save(prop_qSV_fit, file = "data/prop_qSV_fit.Rdata")
# load("data/prop_qSV_fit.Rdata", verbose = TRUE)


#### Tile plots ####
my_breaks <- c(0.05, 0.01, 0.005, 0)

sig_colors <- c(rev(viridis_pal(option = "magma")(3)),NA)
names(sig_colors) <- c("< 0.005","< 0.01", "< 0.05", "> 0.05")


walk2(prop_qSV_fit, names(prop_qSV_fit),function(values, n){
  tile_plot <- values %>%
    ggplot(aes(x = cell_type, y = qSV, fill = p.bonf)) +
    geom_tile(color = "grey") +
    labs(title ="p-values cell-type prop~qSV",
         subtitle = n) +
    geom_text(aes(label = p.bonf.cat, color = p.bonf.cat),size = 3)+
    scale_color_manual(values = sig_colors) +
    scale_fill_viridis(name = "p-value bonf", trans = "log10", option = "magma")
  
  ggsave(filename = paste0("plots/qSV_prop_fit_tileSig-",n,".png"), plot = tile_plot)
})

walk2(prop_qSV_fit, names(prop_qSV_fit),function(values, n){
  tile_plot <- values %>%
    ggplot(aes(x = cell_type, y = qSV, fill = p.bonf)) +
    geom_tile(color = "grey") +
    labs(title ="p-values cell-type prop~qSV",
         subtitle = n) +
    geom_text(aes(label = round(log10(p.bonf),3), color = p.bonf.cat),size = 3)+
    scale_color_manual(values = sig_colors) +
    scale_fill_viridis(name = "p-value bonf", trans = "log10", option = "magma")
  
  ggsave(filename = paste0("plots/qSV_prop_fit_tileVal-",n,".png"), plot = tile_plot)
})


#### Create scatter plots ####
sig_colors <- c(brewer.pal(3, "Set1"),"black")

walk2(est_prop_qsv, names(est_prop_qsv),function(values, n){
  prop_qsv_temp <-  values %>%
    left_join(prop_qSV_fit[[n]] %>% select(cell_type, qSV, p.bonf, p.bonf.cat),
              by = c("cell_type", "qSV"))
  
  scatter_plot <- ggplot(prop_qsv_temp, aes(x = qSV_value, y = prop, color = p.bonf.cat))+
    geom_point(size = .4, alpha = .25) +
    facet_grid(cell_type~qSV, scales = "free")+
    theme_bw(base_size = 10)+
    scale_color_manual(values = sig_colors) +
    theme(legend.text = element_text(size = 15)) + 
    guides(color = guide_legend(override.aes = list(size=5)))
  
  ggsave(filename = paste0("plots/qSV_cellType_scatter-",n,".png"), plot = scatter_plot, width = 26, height = length(unique(prop_qsv_temp$cell_type)))
  message(n)
})


# sgejobs::job_single('deconvo_vs_qSV', create_shell = TRUE, queue= 'bluejay', memory = '5G', command = "Rscript deconvo_vs_qSV.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


# [1] "Reproducibility information:"
# [1] "2021-01-08 18:06:32 EST"
# user  system elapsed 
# 73.682   1.437  76.935 
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value                                      
# version  R version 4.0.3 Patched (2020-11-29 r79529)
# os       CentOS Linux 7 (Core)                      
# system   x86_64, linux-gnu                          
# ui       X11                                        
# language (EN)                                       
# collate  en_US.UTF-8                                
# ctype    en_US.UTF-8                                
# tz       US/Eastern                                 
# date     2021-01-08                                 
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date       lib source                                   
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.0.3)                           
# backports              1.2.0    2020-11-02 [1] CRAN (R 4.0.3)                           
# Biobase              * 2.50.0   2020-10-27 [2] Bioconductor                             
# BiocGenerics         * 0.36.0   2020-10-27 [2] Bioconductor                             
# bitops                 1.0-6    2013-08-17 [2] CRAN (R 4.0.3)                           
# broom                * 0.7.3    2020-12-16 [2] CRAN (R 4.0.3)                           
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.0.3)                           
# cli                    2.2.0    2020-11-20 [1] CRAN (R 4.0.3)                           
# colorspace             2.0-0    2020-11-11 [2] CRAN (R 4.0.3)                           
# crayon                 1.3.4    2017-09-16 [2] CRAN (R 4.0.3)                           
# DBI                    1.1.0    2019-12-15 [2] CRAN (R 4.0.3)                           
# dbplyr                 2.0.0    2020-11-03 [2] CRAN (R 4.0.3)                           
# DelayedArray           0.16.0   2020-10-27 [2] Bioconductor                             
# digest                 0.6.27   2020-10-24 [1] CRAN (R 4.0.3)                           
# dplyr                * 1.0.2    2020-08-18 [1] CRAN (R 4.0.3)                           
# ellipsis               0.3.1    2020-05-15 [2] CRAN (R 4.0.3)                           
# fansi                  0.4.1    2020-01-08 [2] CRAN (R 4.0.3)                           
# farver                 2.0.3    2020-01-16 [2] CRAN (R 4.0.3)                           
# forcats              * 0.5.0    2020-03-01 [2] CRAN (R 4.0.3)                           
# fs                     1.5.0    2020-07-31 [1] CRAN (R 4.0.3)                           
# generics               0.1.0    2020-10-31 [2] CRAN (R 4.0.3)                           
# GenomeInfoDb         * 1.26.2   2020-12-08 [2] Bioconductor                             
# GenomeInfoDbData       1.2.4    2020-11-30 [2] Bioconductor                             
# GenomicRanges        * 1.42.0   2020-10-27 [2] Bioconductor                             
# ggplot2              * 3.3.3    2020-12-30 [2] CRAN (R 4.0.3)                           
# glue                   1.4.2    2020-08-27 [1] CRAN (R 4.0.3)                           
# googledrive            1.0.1    2020-05-05 [1] CRAN (R 4.0.3)                           
# gridExtra              2.3      2017-09-09 [2] CRAN (R 4.0.3)                           
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.0.3)                           
# haven                  2.3.1    2020-06-01 [1] CRAN (R 4.0.3)                           
# here                 * 1.0.0    2020-11-15 [1] CRAN (R 4.0.3)                           
# hms                    0.5.3    2020-01-08 [2] CRAN (R 4.0.3)                           
# httr                   1.4.2    2020-07-20 [1] CRAN (R 4.0.3)                           
# IRanges              * 2.24.0   2020-10-27 [1] Bioconductor                             
# jaffelab             * 0.99.30  2020-11-02 [1] Github (LieberInstitute/jaffelab@42637ff)
# jsonlite               1.7.1    2020-09-07 [1] CRAN (R 4.0.3)                           
# labeling               0.4.2    2020-10-20 [2] CRAN (R 4.0.3)                           
# lattice                0.20-41  2020-04-02 [3] CRAN (R 4.0.3)                           
# lifecycle              0.2.0    2020-03-06 [1] CRAN (R 4.0.3)                           
# limma                  3.46.0   2020-10-27 [2] Bioconductor                             
# lubridate              1.7.9.2  2020-11-13 [1] CRAN (R 4.0.3)                           
# magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.0.3)                           
# Matrix                 1.3-2    2021-01-06 [3] CRAN (R 4.0.3)                           
# MatrixGenerics       * 1.2.0    2020-10-27 [2] Bioconductor                             
# matrixStats          * 0.57.0   2020-09-25 [2] CRAN (R 4.0.3)                           
# modelr                 0.1.8    2020-05-19 [2] CRAN (R 4.0.3)                           
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.0.3)                           
# pillar                 1.4.7    2020-11-20 [1] CRAN (R 4.0.3)                           
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.0.3)                           
# plyr                   1.8.6    2020-03-03 [2] CRAN (R 4.0.3)                           
# ps                     1.4.0    2020-10-07 [1] CRAN (R 4.0.3)                           
# purrr                * 0.3.4    2020-04-17 [1] CRAN (R 4.0.3)                           
# R6                     2.5.0    2020-10-28 [1] CRAN (R 4.0.3)                           
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.0.3)                           
# RColorBrewer         * 1.1-2    2014-12-07 [2] CRAN (R 4.0.3)                           
# Rcpp                   1.0.5    2020-07-06 [1] CRAN (R 4.0.3)                           
# RCurl                  1.98-1.2 2020-04-18 [2] CRAN (R 4.0.3)                           
# readr                * 1.4.0    2020-10-05 [2] CRAN (R 4.0.3)                           
# readxl                 1.3.1    2019-03-13 [2] CRAN (R 4.0.3)                           
# reprex                 0.3.0    2019-05-16 [2] CRAN (R 4.0.3)                           
# reshape2             * 1.4.4    2020-04-09 [2] CRAN (R 4.0.3)                           
# rlang                  0.4.9    2020-11-26 [1] CRAN (R 4.0.3)                           
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.0.3)                           
# rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.0.3)                           
# rvest                  0.3.6    2020-07-25 [2] CRAN (R 4.0.3)                           
# S4Vectors            * 0.28.1   2020-12-09 [2] Bioconductor                             
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.0.3)                           
# segmented              1.3-0    2020-10-27 [1] CRAN (R 4.0.3)                           
# sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.0.3)                           
# stringi                1.5.3    2020-09-09 [1] CRAN (R 4.0.3)                           
# stringr              * 1.4.0    2019-02-10 [2] CRAN (R 4.0.3)                           
# SummarizedExperiment * 1.20.0   2020-10-27 [2] Bioconductor                             
# tibble               * 3.0.4    2020-10-12 [1] CRAN (R 4.0.3)                           
# tidyr                * 1.1.2    2020-08-27 [2] CRAN (R 4.0.3)                           
# tidyselect             1.1.0    2020-05-11 [2] CRAN (R 4.0.3)                           
# tidyverse            * 1.3.0    2019-11-21 [1] CRAN (R 4.0.3)                           
# vctrs                  0.3.5    2020-11-17 [1] CRAN (R 4.0.3)                           
# viridis              * 0.5.1    2018-03-29 [2] CRAN (R 4.0.3)                           
# viridisLite          * 0.3.0    2018-02-01 [2] CRAN (R 4.0.3)                           
# withr                  2.3.0    2020-09-22 [1] CRAN (R 4.0.3)                           
# xml2                   1.3.2    2020-04-23 [2] CRAN (R 4.0.3)                           
# XVector                0.30.0   2020-10-27 [2] Bioconductor                             
# zlibbioc               1.36.0   2020-10-27 [2] Bioconductor                             
# 
# [1] /users/lhuuki/R/4.0.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.0.x/R/4.0.x/lib64/R/library

