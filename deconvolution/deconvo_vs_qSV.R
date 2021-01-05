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
## order celltypes
est_prop_qsv <- map(est_prop_qsv, function(x){
  #alphabetize
  l_alpha <- levels(x$cell_type)[order(levels(x$cell_type))]
  i_inhib <- grep("Inhib",l_alpha)
  i_excit <- grep("Excit",l_alpha)
  l_order <- c(l_alpha[!1:length(l_alpha) %in% c(i_inhib, i_excit)],l_alpha[c(i_excit, i_inhib)])
  levels(x$cell_type) <- l_order
  return(x)
})

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

## Order cell_type factors cell types 
prop_qSV_fit <- map(prop_qSV_fit, function(x){
  #alphabetize
  l_alpha <- levels(x$cell_type)[order(levels(x$cell_type))]
  
  i_inhib <- grep("Inhib",l_alpha)
  i_excit <- grep("Excit",l_alpha)
  
  l_order <- c(l_alpha[!1:length(l_alpha) %in% c(i_inhib, i_excit)],l_alpha[c(i_excit, i_inhib)])
  
  levels(x$cell_type) <- l_order
  return(x)
})
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
