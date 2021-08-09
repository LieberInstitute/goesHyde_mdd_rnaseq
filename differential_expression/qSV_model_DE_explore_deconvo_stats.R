
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

## annotation tables
gene_symbol <- as.data.frame(rowData(rse_gene)) %>% rename(Gene = gencodeID) %>%
  select(Gene, Symbol, ensemblID)

marker_stats2 <- marker_stats %>%
  mutate(marker = rank_ratio <= 25,
         cellType.marker = ifelse(marker, as.character(cellType.target), NA)) %>%
  filter(marker) %>%
  select(ensemblID = gene, cellType.marker)

gene_anno <- left_join(gene_symbol, marker_stats2) %>%
  as_tibble()

## sig colors
sig_colors <- c(RColorBrewer::brewer.pal(3, "Set1"),"dark grey")
names(sig_colors) <- c("sig_Both", "sig_no-deconvo", "sig_deconvo", "None")

n_gene <- nrow(gene_topTables$sacc$no_deconvo)


compare_stats <- function(stat_data, term, other_term ,stat = ".F", p_col = ".adj.P"){
  stat_data <- stat_data %>% select(!contains(other_term))
  
  stat_data_stat <- stat_data %>%
    select(!contains(p_col)) %>%
    pivot_longer(cols = contains(term), names_to = "term", values_to = paste0("deconvo",stat)) %>%
    separate(term, into = c("term","method"), extra = "drop") 
  
  stat_data_p <- stat_data %>% select(!contains(stat)) %>%
    pivot_longer(cols = contains(term), names_to = "term", values_to = "deconvo.adj.P.Val") %>%
    separate(term, into = c("term","method"), extra = "drop")
  
  stat_data2 <- left_join(stat_data_stat, stat_data_p) %>%
    # rename(no_deconvo.adj.P.Val = !!no_deconvo_p, deconvo.adj.P.Val = !!deconvo_p) %>%
    mutate(Signif = case_when(no_deconvo.adj.P.Val < 0.05 & deconvo.adj.P.Val < 0.05 ~"sig_Both",
                              no_deconvo.adj.P.Val < 0.05 ~ "sig_no-deconvo",
                              deconvo.adj.P.Val < 0.05 ~ "sig_deconvo",
                              TRUE ~ "None")) %>%
    left_join(gene_anno)
  
  return(stat_data2)
}

## no-deconvo F vs. deconvo F

outGene_F_stats <- map(gene_topTables, function(dx){
    stats <- do.call("cbind",dx)
    stats <- stats[,grepl("\\.F$|adj.P.Val$", colnames(stats))]
    return(stats)
  })

head(outGene_F_stats$sacc)
dim(outGene_F_stats$sacc)
## GGpairs 
f_plot <- map(outGene_F_stats, ~ggpairs(.x, columns = colnames(.x)[grepl(".F", colnames(.x))],
                                lower = list(continuous = wrap("smooth", alpha = 0.1, size=0.1, color = "blue"))))

walk2(f_plot, names(f_plot), ~ggsave(.x + labs(title = paste("F values:",.y)),
                                     filename = here("differential_expression","plots",paste0("F-stats_ggpairs_",.y,".png"))))

## make long F stats
outGene_F_stats2 <- do.call("rbind",outGene_F_stats)

F_stats <- outGene_F_stats2 %>%
  rownames_to_column("data") %>%
  as_tibble() %>%
  separate(data, into = c("Region","Gene"), extra = "merge") 

F_stats_ilr <- compare_stats(F_stats, term = "ilr", other_term = "prop")
F_stats_ilr %>% count(Region, method, Signif)

F_stats_prop <- compare_stats(F_stats, term = "prop", other_term = "ilr")
F_stats_prop %>% count(Region, method, Signif)

F_stats_deconvo <- list(ilr = F_stats_ilr, prop = F_stats_prop)

F_stats_limits <- map(F_stats_deconvo, ~ .x %>% 
  mutate(signif = cut(deconvo.adj.P.Val, breaks = c(1, 0.05, 0.01, 0),
                      labels = c("<= 0.01", "<= 0.05", "> 0.05"))) %>%
  filter(signif != "> 0.05") %>%
  group_by(signif, Region, method) %>%
  summarize(min_F = min(deconvo.F))
)

F_stats_limits_nd <- F_stats %>% 
  select(Region, Gene, no_deconvo.F, no_deconvo.adj.P.Val)%>% 
  mutate(signif = cut(no_deconvo.adj.P.Val, breaks = c(1, 0.05, 0.01, 0),
                              labels = c("<= 0.01", "<= 0.05", "> 0.05"))) %>%
  filter(signif != "> 0.05") %>%
  group_by(signif, Region) %>%
  summarize(min_F = min(no_deconvo.F))

pwalk(list(f = F_stats_deconvo, l = F_stats_limits, n = names(F_stats_limits)),
      function(f, l, n){
        message(paste0("F-stats_signif-",n,".png"))
        F_stats_prop_scatter <- ggplot(f, aes(x= no_deconvo.F, deconvo.F, color = Signif)) +
          geom_hline(data = l, aes(yintercept = min_F, linetype = signif))+
          geom_vline(data = F_stats_limits_nd, aes(xintercept = min_F, linetype = signif))+
          geom_point(aes(color = Signif), alpha = 0.5, size = 0.5) +
          facet_grid(method ~ Region) +
          scale_color_manual(values = sig_colors)
        
        ggsave(F_stats_prop_scatter, 
               filename = here("differential_expression","plots",paste0("F-stats_signif-",n,".png")),
               width = 10)
        
      })

#### T-statas ####
## compare t stats for each dx coefficent from model with & w/o ILR terms

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
  separate(data, into = c("Region","coef","Gene"), extra = "merge") 


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
  filename = paste0("t-stats_ggpairs_", name,".png")
  ggsave(ggp + labs(title = name), filename = here("differential_expression","plots",filename))
})

#### tables comparing deconvo and no-deconvo t-stats ####

t_stats_ilr <- compare_stats(t_stats, term = "ilr", other_term = "prop", stat = ".t")
t_stats_prop <- compare_stats(t_stats, term = "prop", other_term = "ilr", stat = ".t")

t_stats_ilr %>% group_by(Region, coef, method) %>% count(Signif) %>% filter(Signif != "None")
t_stats_prop %>% group_by(Region, coef, method) %>% count(Signif) %>% filter(Signif != "None")

t_stats_deconvo = list(ilr = t_stats_ilr %>% left_join(deconvo_coef %>% select(-coef_prop)),
                       prop = t_stats_prop %>% left_join(deconvo_coef %>% select(-coef_ilr)))

# save(t_stats, t_stats_ilr, t_stats_prop, file = here("differential_expression","data","t_stats_tables.Rdata"))

t_stats_limits <- map(t_stats_deconvo, ~ .x %>% 
                        mutate(signif = cut(deconvo.adj.P.Val, breaks = c(1, 0.05, 0.01, 0),
                                            labels = c("<= 0.01", "<= 0.05", "> 0.05")),
                               neg = ifelse(deconvo.t < 0, -1, 1)) %>%
                        filter(signif == "<= 0.01") %>%
                        group_by(signif,coef,method, Region, neg) %>%
                        summarize(deconvo.t = min(abs(deconvo.t)) * neg)
)

t_stats_limits_nd <- t_stats %>% 
  select(Region, Gene, no_deconvo.t, no_deconvo.adj.P.Val)%>% 
  mutate(signif = cut(no_deconvo.adj.P.Val, breaks = c(1, 0.05, 0.01, 0),
                      labels = c("<= 0.01", "<= 0.05", "> 0.05"))) %>%
  filter(signif != "> 0.05") %>%
  group_by(signif,Region) %>%
  summarize(min_t = min(abs(no_deconvo.t))) %>%
  mutate(min_t_n = -min_t) %>%
  pivot_longer(cols = c("min_t", "min_t_n"), values_to = "no_deconvo.t")

#### t-stats Plots ####
## simple
pwalk(list(t = t_stats_deconvo, l = t_stats_limits, name = names(t_stats_deconvo)), 
           function(t, l, name){
  
  simple_t_plot <- t %>%
    mutate(alpha = Signif != "None") %>%
    ggplot(aes(x = no_deconvo.t, y = deconvo.t)) +
    geom_hline(data = l, aes(yintercept =  deconvo.t, linetype = signif))+
    geom_vline(data = t_stats_limits_nd %>% filter(signif == "<= 0.01"), 
               aes(xintercept =  no_deconvo.t, linetype = signif))+
    geom_point(aes(color = Signif, alpha = alpha), size = 0.5) +
    facet_grid(method ~ Region + coef) +
    scale_color_manual(values = sig_colors)
  
  ggsave(simple_t_plot ,
         filename = here("differential_expression","plots",paste0("t-stats_signif-",name,".png")),
         width = 15)
  
})

## coef plots
walk2(t_stats_deconvo, names(t_stats_deconvo), function(data, name){
  coef_t_plot <- t_stat_ilr2 %>%
    mutate(alpha = Signif != "None") %>%
    rename(coef_deconvo = starts_with("coef_")) %>%
    ggplot(aes(x = no_deconvo.t, y = deconvo.t)) +
    geom_point(aes(color = Signif, alpha = alpha), size = 0.5) +
    facet_grid(Region + coef ~ method + coef_deconvo)+
    scale_color_manual(values = sig_colors)

  ggsave(coef_t_plot ,
         filename = here("differential_expression","plots",paste0("t-stats_signif_coef-",name,".png")), width = 15)

})

#### prop vs. ilr ####

t_stats_deconvo2 <- map(t_stats_deconvo, ~select(.x, !starts_with("coef_"))) %>%
  do.call("rbind", .)

t_stats_deconvo_term <- t_stats_deconvo2 %>%
  select(Region:deconvo.t) %>%
  pivot_wider(names_from = "term", values_from = "deconvo.t",
              names_prefix = "t.") %>%
  left_join(t_stats_deconvo2 %>%
              select(Region:method, deconvo.adj.P.Val) %>%
              pivot_wider(names_from = "term", values_from = "deconvo.adj.P.Val",
                          names_prefix = "adj.P.Val.")) %>%
  mutate(Signif = case_when(adj.P.Val.ilr < 0.05 & adj.P.Val.prop < 0.05 ~"sig_Both",
                            adj.P.Val.ilr < 0.05 ~ "sig_ilr",
                            adj.P.Val.prop < 0.05 ~ "sig_prop",
                            TRUE ~ "None")) %>%
  left_join(gene_anno)

## sig colors
names(sig_colors) <- c("sig_Both", "sig_ilr", "sig_prop", "None")

simple_t_plot_term <- t_stats_deconvo_term  %>%
  mutate(alpha = Signif != "None") %>%
  ggplot(aes(x = t.prop, y = t.ilr)) +
  geom_point(aes(color = Signif, alpha = alpha), size = 0.5) +
  facet_grid(method ~ Region + coef) +
  scale_color_manual(values = sig_colors)

ggsave(simple_t_plot_term ,
       filename = here("differential_expression","plots","t-stats_signif-terms.png"),
       width = 15)

### by method
t_stats_deconvo_method <- t_stats_deconvo2 %>%
  select(Region:deconvo.t) %>%
  pivot_wider(names_from = "method", values_from = "deconvo.t",
              names_prefix = "t.") %>%
  left_join(t_stats_deconvo2 %>%
              select(Region:method, deconvo.adj.P.Val) %>%
              pivot_wider(names_from = "method", values_from = "deconvo.adj.P.Val",
                          names_prefix = "adj.P.Val.")) %>%
  mutate(Signif = case_when(adj.P.Val.music < 0.05 & adj.P.Val.bisque < 0.05 ~"sig_Both",
                            adj.P.Val.music < 0.05 ~ "sig_music",
                            adj.P.Val.bisque < 0.05 ~ "sig_bisque",
                            TRUE ~ "None")) %>%
  left_join(gene_anno)

t_stats_deconvo_method %>% count(Signif)

names(sig_colors) <- c("sig_Both", "sig_music", "sig_bisque", "None")

simple_t_plot_method <- t_stats_deconvo_method  %>%
  mutate(alpha = Signif != "None") %>%
  ggplot(aes(x = t.bisque, y = t.music)) +
  geom_point(aes(color = Signif, alpha = alpha), size = 0.5) +
  facet_grid(term ~ Region + coef) +
  scale_color_manual(values = sig_colors)

ggsave(simple_t_plot_method ,
       filename = here("differential_expression","plots","t-stats_signif-method.png"),
       width = 15)

#sgejobs::job_single('qSV_model_DE_explore', create_shell = TRUE, memory = '80G', command = "Rscript qSV_model_DE_explore.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
