library(ggplot2)
library(ggrepel)
library(dplyr)
library(data.table)
library(plotly)
library(htmlwidgets)
library(sessioninfo)
library(SummarizedExperiment)
library(ggpubr)
library(tools)
library(xlsx)

data.table::setDTthreads(threads = 1)

# Sourcing Data/Inst. Vars. ####
load("rda/twas_exp_ranges.Rdata")
# load("twas_exp_ranges.Rdata")

dir.create(file.path("analysis/plots"),
           showWarnings = FALSE,
           recursive = TRUE)
dir.create(file.path("analysis/tables"),
           showWarnings = FALSE,
           recursive = TRUE)

# Filter N/A Z scores
twas_z <- twas_exp_fin %>% filter(!is.na(TWAS.Z))

twas_z_sACC <- twas_z[twas_z$region == "sACC", ]

twas_z_amyg <- twas_z[twas_z$region == "Amygdala", ]

don <- list()

axisdf <- list()

don_key <- list()

p <- list()

intctv_plot <- list()

fin_plot <- list()

# Preprocessing Data ####
for (i in 1:2) {
    if (i == 1) {
        twas_var <- twas_z_amyg
    } else {
        twas_var <- twas_z_sACC
    }
    
    don[[i]] <-
        twas_var %>%
        # Compute chromosome size
        group_by(CHR) %>%
        summarise(chr_len = max(end)) %>%
        
        # Calculate cumulative position of each chromosome
        mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
        select(-chr_len) %>%
        
        # Add this info to the initial dataset
        left_join(twas_var,
                  .,
                  by = c("CHR" = "CHR")) %>%
        
        # Add a cumulative position of each SNP
        arrange(CHR, twas_mean_dist) %>%
        mutate(BPcum = twas_mean_dist + tot)
    
    
    axisdf[[i]] = don[[i]] %>% group_by(CHR) %>% summarise(center = (max(BPcum) + min(BPcum)) / 2)
    
    # Prepare text description for each SNP:
    don[[i]]$text <-
        paste0(
            "Gene Symbol: ",
            don[[i]]$genesymbol,
            "\nENSEMBL Gene ID: ",
            don[[i]]$geneid,
            "\nBrain Subregion: ",
            don[[i]]$region,
            "\nChromosome: ",
            don[[i]]$CHR,
            "\nStart Position: ",
            don[[i]]$start,
            "\nEnd Position: ",
            don[[i]]$end,
            "\nZ score: ",
            don[[i]]$TWAS.Z %>% round(2)
        )
    
    don_key[[i]] <-
        highlight_key(don[[i]], ~ genesymbol, group = "Gene Symbol")
}

# TWAS Z Manhattan Plot ####
pdf(file = here::here("twas_both", "analysis", "plots", "goesHyde_TWAS_ManhattanPlot.pdf"))
# storing ggplot as an object3

sig <- qnorm(1 - 0.025 / table(twas_exp_fin$region))
for (i in 1:2) {
    # Bonferroni Correction
    sig_bonf <- sig[[i]]
    
    p[[i]] <-
        ggplot(don_key[[i]], aes(x = BPcum, y = TWAS.Z, text = text)) +
        
        ggtitle(paste0("Gene Windows of ", ifelse(i == 1, "Amygdala", "sACC") , " TWAS")) +
        # Show all points
        geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
        scale_color_manual(values = rep(c("#861657", "#D56AA0"), 22)) +
        geom_hline(
            yintercept = c(sig_bonf,-sig_bonf),
            color = "grey40",
            linetype = "dashed"
        ) +
        
        # custom X axis:
        scale_x_continuous(labels = axisdf[[i]]$CHR, breaks = axisdf[[i]]$center) +
        scale_y_continuous(expand = c(0, 0)) +     # remove space between plot area and x axis
        
        # Custom the theme:
        theme_bw() +
        theme(
            legend.position = "none",
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()
        )
    
    print(p[[i]])
}
dev.off()

# Z scores threshold
twas_z_amyg_threshold <-
    rbind(twas_z_amyg[TWAS.Z > sig[[1]], ], twas_z_amyg[TWAS.Z < -sig[[1]], ])
twas_z_sACC_threshold <-
    rbind(twas_z_sACC[TWAS.Z > sig[[2]], ], twas_z_sACC[TWAS.Z < -sig[[2]], ])

# Interactive TWAS Z Manhattan Plots ####
for (i in 1:2) {
    ##### Plotly
    intctv_plot[[i]] <- ggplotly(p[[i]], tooltip = "text")
    
    fin_plot[[i]] <- highlight(
        intctv_plot[[i]],
        on = "plotly_click",
        off = "plotly_doubleclick",
        color = "#60D394",
        selectize = TRUE
    )
    
    saveWidget(fin_plot[[i]],
               file.path(paste0(
                   # "analysis/plots/",
                   "goesHyde_TWAS_",
                   ifelse(i == 1, "Amygdala", "sACC"),
                   "_ManhattanPlotly.html"
               )))
}

# Plotly cannot save directly to relative path for whatever reason
system("mv *_ManhattanPlotly.html analysis/plots/")

# Scatter plots ####

pdf(
    'analysis/plots/goesHyde_TWAS_ScatterPlots_hsqFail.pdf',
    useDingbats = FALSE,
    width = 10,
    height = 10
)

twas_z_select <-
    select(twas_z, geneid, genesymbol, TWAS.Z, TWAS.P, region) %>%
    as.data.table()

# Render Z scores and P values horizontally by region
twas_z_wide <-
    dcast(twas_z_select,
          geneid + genesymbol ~ region,
          value.var = c("TWAS.Z", "TWAS.P"))

# FDR calculation per subregion
twas_z_wide$Amygdala.fdr.p <-
    p.adjust(twas_z_wide$TWAS.P_Amygdala, 'fdr')
twas_z_wide$sACC.fdr.p <- p.adjust(twas_z_wide$TWAS.P_sACC, 'fdr')

# Indicate in both
# twas_z_wide$in_both <-
#     ifelse(!is.na(twas_z_wide$TWAS.Z_Amygdala &
#                       twas_z_wide$TWAS.Z_sACC),
#            TRUE,
#            FALSE)

# FDR cutoffs
twas_z_wide$FDR.5perc <- 'None'
twas_z_wide$FDR.5perc[twas_z_wide$Amygdala.fdr.p < 0.05] <-
    'Amygdala'
twas_z_wide$FDR.5perc[twas_z_wide$sACC.fdr.p < 0.05] <- 'sACC'
twas_z_wide$FDR.5perc[twas_z_wide$Amygdala.fdr.p < 0.05 &
                          twas_z_wide$sACC.fdr.p < 0.05] <- 'Both'

# Remove NAs
twas_z_wide[is.na(twas_z_wide)] <- 0
# twas_z_wide <- twas_z_wide[twas_z_wide$in_both, ]

twas_z_wide$FDR.5perc <-
    factor(twas_z_wide$FDR.5perc,
           levels = c('None', 'Amygdala', 'sACC', 'Both'))

ggplot(twas_z_wide,
       aes(x = TWAS.Z_Amygdala,
           y = TWAS.Z_sACC,
           color = FDR.5perc)) +
    xlab("Amygdala Z-Score") +
    ylab("sACC Z-Score") +
    labs(color = "FDR < 5%") +
    geom_point() +
    coord_fixed() +
    theme_bw(base_size = 20) +
    scale_color_manual(values = c('grey80', 'dark orange', 'skyblue3', 'purple')) # you can define names

dev.off()

## Z-Score Correlation Test ####

both_genes_sACC <-
    twas_z_sACC[twas_z_sACC$geneid %in% twas_z_amyg$geneid]
both_genes_amyg <-
    twas_z_amyg[twas_z_amyg$geneid %in% twas_z_sACC$geneid]

z_score_cor <-
    cor.test(both_genes_sACC$TWAS.Z, both_genes_amyg$TWAS.Z, method = "pearson")
z_score_cor

both_z_scores <- both_genes_sACC %>%
    select(TWAS.Z) %>%
    mutate(
        sACC_TWAS_Z_both = TWAS.Z,
        amyg_TWAS_Z_both = both_genes_amyg$TWAS.Z,
        .keep = "unused"
    )

pdf(
    'analysis/plots/goesHyde_TWAS_Z_Correlation.pdf',
    useDingbats = FALSE,
    width = 10,
    height = 10
)


ggscatter(
    both_z_scores,
    x = "amyg_TWAS_Z_both",
    y = "sACC_TWAS_Z_both",
    add = "reg.line",
    add.params = list(color = "red", fill = "orangered1"),
    conf.int = TRUE,
    cor.coef = TRUE,
    cor.method = "pearson",
    alpha = 0.5,
    xlab = toTitleCase("TWAS Z-Scores for Amygdala genes"),
    ylab = toTitleCase("TWAS Z-Scores for sACC genes")
) +
    grids(size = 1) + bgcolor("gray90") + border("gray90")

dev.off()

# Differential Expression ####
# Skipping diff exp for now because I don't have the corresponding file
if (FALSE) {
    load(
        "/dcl01/lieber/ajaffe/lab/zandiHyde_bipolar_rnaseq/case_control/bipolarControl_deStats_byRegion_qSVAjoint_withAnnotation.rda",
        verbose = TRUE
    )
    
    statOutGene <- statOut[statOut$Type == "Gene", ] %>%
        as.data.table(keep.rownames = "geneid") %>%
        select(geneid, t_Amyg, t_sACC)
    
    
    merged_t <- merge(statOutGene, twas_z_wide, by = "geneid")
    
    pdf(
        'analysis/plots/MDD_TWAS_t-stat.pdf',
        useDingbats = FALSE,
        width = 10,
        height = 10
    )
    
    amyg_rho <-
        format(
            cor(merged_t$TWAS.Z_Amygdala, merged_t$t_Amyg),
            digits = 3,
            scientific = TRUE
        )
    sACC_rho <-
        format(
            cor(merged_t$TWAS.Z_sACC, merged_t$t_sACC),
            digits = 3,
            scientific = TRUE
        )
    
    ggplot(merged_t,
           aes(x = TWAS.Z_Amygdala,
               y = t_Amyg)) + geom_point() + labs(
                   title = toTitleCase("TWAS vs MDD differential expression in Amygdala"),
                   x = "TWAS Z-Score",
                   y = "MDD vs Control t-Statistic"
               ) + annotate(
                   "text",
                   x = -5.5,
                   y = 6,
                   label = paste0("rho == ", formatC(amyg_rho, format = "e")),
                   parse = TRUE
               ) + scale_y_continuous(breaks = c(-6,-3, 0, 3, 6)) + xlim(-6, 6) +
        theme_bw(base_size = 20)
    
    ggplot(merged_t,
           aes(x = TWAS.Z_sACC,
               y = t_sACC)) + geom_point() + labs(
                   title = toTitleCase("TWAS vs MDD differential expression in sACC"),
                   x = "TWAS Z-Score",
                   y = "MDD vs Control t-Statistic"
               ) + annotate(
                   "text",
                   x = -5.5,
                   y = 6,
                   label = paste0("rho == ", formatC(sACC_rho, format = "e")),
                   parse = TRUE
               ) + scale_y_continuous(breaks = c(-6,-3, 0, 3, 6)) +
        scale_x_continuous(breaks = waiver()) +
        theme_bw(base_size = 20)
    
    dev.off()
}

# XLSX Output ####
write.xlsx2(
    x = twas_z_amyg_threshold,
    file = "analysis/tables/goesHyde_Amyg_sACC_FinalOutputTable.xlsx",
    sheetName = "Significant TWAS Z Scores in Amygdala",
    col.names = TRUE,
    row.names = FALSE,
    append = FALSE
)

write.xlsx2(
    x = twas_z_sACC_threshold,
    file = "analysis/tables/goesHyde_Amyg_sACC_FinalOutputTable.xlsx",
    sheetName = "Significant TWAS Z Scores in sACC",
    col.names = TRUE,
    row.names = FALSE,
    append = TRUE
)

write.xlsx2(
    x = twas_z_wide,
    file = "analysis/tables/goesHyde_Amyg_sACC_FinalOutputTable.xlsx",
    sheetName = "TWAS Z Scatterplot with FDR and P-Values for Both Regions",
    col.names = TRUE,
    row.names = FALSE,
    append = TRUE
)

# write.xlsx2(x = merged_t, file = "analysis/tables/MDD_Amyg_sACC_FinalOutputTable.xlsx", sheetName = "TWAS vs MDD Differential Expression in Both Regions", col.names = TRUE, row.names = FALSE, append = TRUE)

# for enrichment test
save.image("rda/generate_plots_data.RData")

## Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# EXECUTED INTERACTIVELY
# > print("Reproducibility information:")
# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2021-09-30 12:28:57 EDT"
# > proc.time()
# user  system elapsed 
# 59.658   5.185 742.482 
# > options(width = 120)
# > session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value                                      
# version  R version 4.1.0 Patched (2021-05-18 r80330)
# os       CentOS Linux 7 (Core)                      
# system   x86_64, linux-gnu                          
# ui       X11                                        
# language (EN)                                       
# collate  en_US.UTF-8                                
# ctype    en_US.UTF-8                                
# tz       US/Eastern                                 
# date     2021-09-30                                 
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date       lib source        
# abind                  1.4-5    2016-07-21 [2] CRAN (R 4.1.0)
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# backports              1.2.1    2020-12-09 [2] CRAN (R 4.1.0)
# Biobase              * 2.52.0   2021-05-19 [2] Bioconductor  
# BiocGenerics         * 0.38.0   2021-05-19 [2] Bioconductor  
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# broom                  0.7.9    2021-07-27 [2] CRAN (R 4.1.0)
# car                    3.0-11   2021-06-27 [2] CRAN (R 4.1.0)
# carData                3.0-4    2020-05-22 [2] CRAN (R 4.1.0)
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.1.0)
# cli                    3.0.1    2021-07-17 [2] CRAN (R 4.1.0)
# colorspace             2.0-2    2021-06-24 [2] CRAN (R 4.1.0)
# crayon                 1.4.1    2021-02-08 [2] CRAN (R 4.1.0)
# crosstalk              1.1.1    2021-01-12 [2] CRAN (R 4.1.0)
# curl                   4.3.2    2021-06-23 [2] CRAN (R 4.1.0)
# data.table           * 1.14.0   2021-02-21 [2] CRAN (R 4.1.0)
# DBI                    1.1.1    2021-01-15 [2] CRAN (R 4.1.0)
# DelayedArray           0.18.0   2021-05-19 [2] Bioconductor  
# digest                 0.6.27   2020-10-24 [2] CRAN (R 4.1.0)
# dplyr                * 1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  0.5.0    2021-05-25 [2] CRAN (R 4.1.0)
# farver                 2.1.0    2021-02-28 [2] CRAN (R 4.1.0)
# fastmap                1.1.0    2021-01-25 [2] CRAN (R 4.1.0)
# forcats                0.5.1    2021-01-27 [2] CRAN (R 4.1.0)
# foreign                0.8-81   2020-12-22 [3] CRAN (R 4.1.0)
# generics               0.1.0    2020-10-31 [2] CRAN (R 4.1.0)
# GenomeInfoDb         * 1.28.1   2021-07-01 [2] Bioconductor  
# GenomeInfoDbData       1.2.6    2021-05-11 [2] Bioconductor  
# GenomicRanges        * 1.44.0   2021-05-19 [2] Bioconductor  
# ggplot2              * 3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
# ggpubr               * 0.4.0    2020-06-27 [1] CRAN (R 4.1.0)
# ggrepel              * 0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
# ggsignif               0.6.3    2021-09-09 [1] CRAN (R 4.1.0)
# glue                   1.4.2    2020-08-27 [2] CRAN (R 4.1.0)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# haven                  2.4.3    2021-08-04 [2] CRAN (R 4.1.0)
# here                   1.0.1    2020-12-13 [1] CRAN (R 4.1.0)
# hms                    1.1.0    2021-05-17 [2] CRAN (R 4.1.0)
# htmltools              0.5.1.1  2021-01-22 [2] CRAN (R 4.1.0)
# htmlwidgets          * 1.5.3    2020-12-10 [2] CRAN (R 4.1.0)
# httpuv                 1.6.1    2021-05-07 [2] CRAN (R 4.1.0)
# httr                   1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
# IRanges              * 2.26.0   2021-05-19 [2] Bioconductor  
# jsonlite               1.7.2    2020-12-09 [2] CRAN (R 4.1.0)
# labeling               0.4.2    2020-10-20 [2] CRAN (R 4.1.0)
# later                  1.2.0    2021-04-23 [2] CRAN (R 4.1.0)
# lattice                0.20-44  2021-05-02 [3] CRAN (R 4.1.0)
# lazyeval               0.2.2    2019-03-15 [2] CRAN (R 4.1.0)
# lifecycle              1.0.0    2021-02-15 [2] CRAN (R 4.1.0)
# magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
# Matrix                 1.3-4    2021-06-01 [3] CRAN (R 4.1.0)
# MatrixGenerics       * 1.4.2    2021-08-08 [2] Bioconductor  
# matrixStats          * 0.60.0   2021-07-26 [2] CRAN (R 4.1.0)
# mgcv                   1.8-36   2021-06-01 [3] CRAN (R 4.1.0)
# mime                   0.11     2021-06-23 [2] CRAN (R 4.1.0)
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# nlme                   3.1-152  2021-02-04 [3] CRAN (R 4.1.0)
# openxlsx               4.2.4    2021-06-16 [2] CRAN (R 4.1.0)
# pillar                 1.6.2    2021-07-29 [2] CRAN (R 4.1.0)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# plotly               * 4.9.4.1  2021-06-18 [2] CRAN (R 4.1.0)
# promises               1.2.0.1  2021-02-11 [2] CRAN (R 4.1.0)
# purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R6                     2.5.0    2020-10-28 [2] CRAN (R 4.1.0)
# Rcpp                   1.0.7    2021-07-07 [2] CRAN (R 4.1.0)
# RCurl                  1.98-1.3 2021-03-16 [2] CRAN (R 4.1.0)
# readxl                 1.3.1    2019-03-13 [2] CRAN (R 4.1.0)
# rio                    0.5.27   2021-06-21 [2] CRAN (R 4.1.0)
# rJava                  1.0-4    2021-04-29 [2] CRAN (R 4.1.0)
# rlang                  0.4.11   2021-04-30 [2] CRAN (R 4.1.0)
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
# rstatix                0.7.0    2021-02-13 [1] CRAN (R 4.1.0)
# S4Vectors            * 0.30.0   2021-05-19 [2] Bioconductor  
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
# sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.1.0)
# shiny                  1.6.0    2021-01-25 [2] CRAN (R 4.1.0)
# stringi                1.7.3    2021-07-16 [2] CRAN (R 4.1.0)
# SummarizedExperiment * 1.22.0   2021-05-19 [2] Bioconductor  
# tibble                 3.1.3    2021-07-23 [2] CRAN (R 4.1.0)
# tidyr                  1.1.3    2021-03-03 [2] CRAN (R 4.1.0)
# tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
# viridisLite            0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
# withr                  2.4.2    2021-04-18 [2] CRAN (R 4.1.0)
# xlsx                 * 0.6.5    2020-11-10 [2] CRAN (R 4.1.0)
# xlsxjars               0.6.1    2014-08-22 [2] CRAN (R 4.1.0)
# xtable                 1.8-4    2019-04-21 [2] CRAN (R 4.1.0)
# XVector                0.32.0   2021-05-19 [2] Bioconductor  
# yaml                   2.2.1    2020-02-01 [2] CRAN (R 4.1.0)
# zip                    2.2.0    2021-05-31 [2] CRAN (R 4.1.0)
# zlibbioc               1.38.0   2021-05-19 [2] Bioconductor  
# 
# [1] /users/aseyedia/R/4.1
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/library
