
library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(devtools)
library(sessioninfo)
library(tidyverse)
library(GGally)
library(viridis)
library(plotly)
library(here)

#### Load Data ####

load(here("deconvolution","data","cell_colors.Rdata"), verbose= TRUE)
load(here("differential_expression","data","qSVA_MDD_gene_DEresults.rda"), verbose = TRUE)
load(here("differential_expression","data","deconvo_coef_explore.Rdata"), verbose = TRUE)
load(here("deconvolution","data","marker_stats.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)

gene_symbol <- as.data.frame(rowData(rse_gene)) %>% rename(Gene = gencodeID) %>%
  select(Gene, Symbol)

n_gene <- nrow(gene_topTables$sacc$no_deconvo)

f_table <- map(gene_topTables, function(tt){
  tt_all_col <- do.call("cbind",tt)
  tt_F <- tt_all_col[,grepl(".F",colnames(tt_all_col))]
  return(tt_F)
  })

f_plot <- map(f_table, ~ggpairs(.x,lower = list(continuous = wrap("smooth", alpha = 0.1, size=0.1, color = "blue"))))

walk2(f_plot, names(f_plot), ~ggsave(.x + labs(title = paste("F values:",.y)),
                                     filename = here("differential_expression","plots",paste0("F_plot-",.y,".png"))))

## compare t stats for each dx coefficent from model with & w/o ILR terms
stat_cols <- c("t","adj.P.Val")
outGene_t_stats <- map(outGene_single_coef, function(region){
  dx_stats <- map(region, function(dx){
    stats <- do.call("cbind",dx)
    stats <- stats[,grepl("\\.t$|adj.P.Val$", colnames(stats))]
    return(stats)
  })
  return(dx_stats)
})

head(outGene_t_stats$sacc$ctrl)

outGene_t_stats2 <- map(outGene_t_stats, ~do.call("rbind",.x))
outGene_t_stats3 <- do.call("rbind", outGene_t_stats2)

deconvo_coef <- do.call("rbind", combo_coef_explore) %>%
  add_column(Region = rep(names(combo_coef_explore), each = nrow(combo_coef_explore[[1]]))) %>%
  rename(Gene = gencodeID)

deconvo_coef %>% count(Region, method)

t_stats <- outGene_t_stats3 %>%
  rownames_to_column("data") %>%
  as_tibble() %>%
  separate(data, into = c("Region","coef","Gene"), extra = "merge") %>%
  left_join(combo_coef_explore2)


## ggpairs for t-stats
t_stats_grouped_ggpairs <- t_stats %>%
  group_by(Region, coef) %>%
  group_map(~ggpairs(.x, c("no_deconvo.t", "prop_music.t","prop_bisque.t"),
                     lower = list(continuous = wrap("smooth", alpha = 0.1, size=0.1, color = "blue"))))

names(t_stats_grouped_ggpairs) <- t_stats %>%
  group_by(Region, coef) %>%
  count() %>%
  mutate(label = paste0(Region,"_",coef)) %>%
  pull(label)

walk2(t_stats_grouped_ggpairs, names(t_stats_grouped_ggpairs), function(ggp, name){
  message(name)
  filename = paste0("ggpair_t-stats_", name,".png")
  ggsave(ggp + labs(title = name), filename = here("differential_expression","plots",filename))
})

#### tables comparing deconvo and no-deconvo t-stats ####
marker_stats2 <- marker_stats %>%
  mutate(marker = rank_ratio <= 25,
         cellType.marker = ifelse(marker, as.character(cellType.target), NA)) %>%
  filter(marker) %>%
  select(ensemblID = gene, cellType.marker)

t_stats_ilr <-  t_stats %>%
  select(!contains("prop"))

t_stats_ilr2 <-  t_stats_ilr %>% select(!contains("adj.P")) %>%
  pivot_longer(cols = contains("ilr"), names_to = "term", values_to = "deconvo.t") %>%
  separate(term, into = c("term","method"), extra = "drop") %>%
  left_join(t_stats_ilr %>% select(!contains(".t")) %>%
              pivot_longer(cols = contains("ilr"), names_to = "term", values_to = "deconvo.adj.P.Val") %>%
              separate(term, into = c("term","method"), extra = "drop"))

t_stats_ilr3 <- t_stats_ilr2 %>%
    mutate(Signif = case_when(no_deconvo.adj.P.Val < 0.05 & deconvo.adj.P.Val < 0.05 ~"sig_Both",
                              no_deconvo.adj.P.Val < 0.05 ~ "sig_no-deconvo",
                              deconvo.adj.P.Val < 0.05 ~ "sig_deconvo-ilr",
                              TRUE ~ "None"),
           ensemblID = ss(Gene, "\\.")) %>%
  left_join(deconvo_coef %>% select(- coef_prop)) %>%
  left_join(marker_stats2) %>%
  left_join(gene_symbol)

t_stats_ilr3 %>% group_by(Region, coef, method) %>% count(Signif, cellType.marker) %>% filter(Signif != "None")

## prop
t_stats_prop <- t_stats %>%
  select(!contains("ilr"))

t_stats_prop2 <-  t_stats_prop %>%
  select(!contains("ilr")) %>%
  select(!contains("adj.P")) %>%
  pivot_longer(cols = contains("prop"), names_to = "term", values_to = "deconvo.t") %>%
  separate(term, into = c("term","method"), extra = "drop") %>%
  left_join(t_stats_prop %>% select(!contains(".t")) %>%
              pivot_longer(cols = contains("prop"), names_to = "term", values_to = "deconvo.adj.P.Val") %>%
              separate(term, into = c("term","method"), extra = "drop"))

t_stats_prop3 <- t_stats_prop2 %>%
  mutate(Signif = case_when(no_deconvo.adj.P.Val < 0.05 & deconvo.adj.P.Val < 0.05 ~"sig_Both",
                            no_deconvo.adj.P.Val < 0.05 ~ "sig_no-deconvo",
                            deconvo.adj.P.Val < 0.05 ~ "sig_deconvo-prop",
                            TRUE ~ "None"),
         ensemblID = ss(Gene, "\\.")) %>%
  left_join(deconvo_coef %>% select(- coef_ilr))%>%
  left_join(marker_stats2)%>%
  left_join(gene_symbol)

t_stats_prop3 %>% group_by(Region, coef, method) %>% count(Signif, cellType.marker) %>% filter(Signif != "None")

save(t_stats, t_stats_ilr3, t_stats_prop3, file = here("differential_expression","data","t_stats_tables.Rdata"))

#### t-stats Plots ####
## sig colors
sig_colors <- c(RColorBrewer::brewer.pal(3, "Set1"),"dark grey")
names(sig_colors) <- c("sig_Both", "sig_no-deconvo", "sig_deconvo-prop", "None")

## ilr plots
ilr_plot <- t_stats_ilr3 %>%
  mutate(alpha = Signif != "None") %>%
  ggplot(aes(x = no_deconvo.t, y = deconvo.t)) +
  geom_point(aes(color = Signif, alpha = alpha), size = 0.5) +
  facet_grid(Region + coef ~ method + coef_ilr)+
  scale_color_manual(values = sig_colors)

ggsave(ilr_plot ,
       filename = here("differential_expression","plots","ilr_plots.png"), width = 15)

## Prop plots
simple_prop_plot <- t_stats_prop3 %>%
  mutate(alpha = Signif != "None") %>%
  ggplot(aes(x = no_deconvo.t, y = deconvo.t)) +
  geom_point(aes(color = Signif,alpha = alpha), size = 0.5) +
  facet_grid(method ~ Region + coef) +
  scale_color_manual(values = sig_colors)

ggsave(simple_prop_plot ,
       filename = here("differential_expression","plots","prop_plots_simple.png"), width = 15)

simple_prop_plot <- t_stats_prop3 %>%
  mutate(alpha = Signif != "None") %>%
  ggplot(aes(x = no_deconvo.t, y = deconvo.t)) +
  geom_point(aes(color = cellType.marker), size = 0.5) +
  facet_grid(method ~ Region + coef) +
  scale_color_manual(values = cell_colors)

ggsave(simple_prop_plot ,
       filename = here("differential_expression","plots","prop_plots_simple-marker.png"), width = 15)


prop_plot <- t_stats_prop3 %>%
  mutate(alpha = Signif != "None") %>%
  ggplot(aes(x = no_deconvo.t, y = deconvo.t)) +
  facet_grid(Region + coef ~ method + coef_prop) +
  geom_point(aes(color = Signif, alpha = alpha), size = 0.5) +
  scale_color_manual(values = sig_colors)

ggsave(prop_plot ,
       filename = here("differential_expression","plots","prop_plots.png"), width = 15)

#sgejobs::job_single('qSV_model_DE_explore', create_shell = TRUE, memory = '80G', command = "Rscript qSV_model_DE_explore.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
