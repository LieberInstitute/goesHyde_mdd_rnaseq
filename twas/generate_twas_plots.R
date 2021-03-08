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

twas_z_sACC <- twas_z[twas_z$region == "sACC",]

twas_z_amyg <- twas_z[twas_z$region == "Amygdala",]

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
    } else{
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
pdf(file = "analysis/plots/MDD_TWAS_ManhattanPlot.pdf")
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
            yintercept = c(sig_bonf, -sig_bonf),
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
    rbind(twas_z_amyg[TWAS.Z > sig[[1]],], twas_z_amyg[TWAS.Z < -sig[[1]],])
twas_z_sACC_threshold <-
    rbind(twas_z_sACC[TWAS.Z > sig[[2]],], twas_z_sACC[TWAS.Z < -sig[[2]],])

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
                   "MDD_TWAS_",
                   ifelse(i == 1, "Amygdala", "sACC"),
                   "_ManhattanPlotly.html"
               )))
}

# Plotly cannot save directly to relative path for whatever reason
system("mv *_ManhattanPlotly.html analysis/plots/")

# Scatter plots ####

pdf(
    'analysis/plots/MDD_TWAS_ScatterPlots.pdf',
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
twas_z_wide$in_both <-
    ifelse(!is.na(twas_z_wide$TWAS.Z_Amygdala &
                      twas_z_wide$TWAS.Z_sACC),
           TRUE,
           FALSE)

# FDR cutoffs
twas_z_wide$FDR.5perc <- 'None'
twas_z_wide$FDR.5perc[twas_z_wide$Amygdala.fdr.p < 0.05] <-
    'Amygdala'
twas_z_wide$FDR.5perc[twas_z_wide$sACC.fdr.p < 0.05] <- 'sACC'
twas_z_wide$FDR.5perc[twas_z_wide$Amygdala.fdr.p < 0.05 &
                          twas_z_wide$sACC.fdr.p < 0.05] <- 'Both'

# Remove NAs
twas_z_wide[is.na(twas_z_wide)] <- 0

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
    'analysis/plots/MDD_TWAS_Z_Correlation.pdf',
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
    
    statOutGene <- statOut[statOut$Type == "Gene",] %>%
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
               ) + scale_y_continuous(breaks = c(-6, -3, 0, 3, 6)) + xlim(-6, 6) +
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
               ) + scale_y_continuous(breaks = c(-6, -3, 0, 3, 6)) +
        scale_x_continuous(breaks = waiver()) +
        theme_bw(base_size = 20)
    
    dev.off()
    
    # XLSX Output ####
    
    # xlsx does not work with conda_R/4.0, needs 4.0.x
    save(
        twas_z_amyg_threshold,
        twas_z_sACC_threshold,
        twas_z_wide,
        merged_t,
        "rda/xlsx_output.RData"
    )
}
# for enrichment test
save.image("rda/generate_plots_data.RData")

## Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
