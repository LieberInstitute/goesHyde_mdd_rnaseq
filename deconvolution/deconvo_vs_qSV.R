library(SummarizedExperiment)
library(jaffelab)
library(here)
library(reshape2)
library(purrr)
library(tidyverse)
library(broom)
library(sessioninfo)

#### Load Data ####
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
pd <- as.data.frame(colData(rse_gene))

## Bind with qSV table
load(here("differential_expression" ,"qSV_mat.Rdata"), verbose = TRUE)
all(rownames(pd)==rownames(qSV_mat))
pd <- cbind(pd, qSV_mat)

## split by region
pd_sacc <- pd[pd$BrainRegion == "sACC",]
pd_amyg <- pd[pd$BrainRegion == "Amygdala",]


#### check against qSV ####
qsv_names <- colnames(qSV_mat)

qsv_long_sacc <- pd_sacc %>% 
  rownames_to_column("sample") %>%
  select(all_of(c("sample", qsv_names))) %>%
  melt(id.vars= c("sample")) %>%
  rename(qSV_name = variable, qSV = value) 

qsv_long_amyg <- pd_amyg %>% 
  rownames_to_column("sample") %>%
  select(all_of(c("sample", qsv_names))) %>%
  melt(id.vars= c("sample")) %>%
  rename(qSV_name = variable, qSV = value) 

prop_qsv_sacc <- prop_sacc %>%
  left_join(qsv_long_sacc, by = "sample")

prop_qsv_amyg <- prop_amyg %>%
  left_join(qsv_long_amyg, by = "sample")

#### Calculate p-values ####
##sACC
cross_sacc <- cross2(cells_sacc, qsv_names)

cross_df_sacc <- as.data.frame(do.call("rbind", cross_sacc))
colnames(cross_df_sacc) <- c("cell_type","qSV_name")
cross_df_sacc$cell_type <- factor(cross_df_sacc$cell_type, levels = unique(cross_df_sacc$cell_type))
cross_df_sacc$qSV_name <- factor(cross_df_sacc$qSV_name, levels = unique(cross_df_sacc$qSV_name))

## All genes
qsv_lm_sacc <- bind_rows(map(cross_sacc, function(x){
  cell = x[[1]]
  qsv = x[[2]]
  prop_qsv <- prop_qsv_sacc %>%
    filter(cell_type == cell & qSV_name == qsv) 
  
  tidy(lm(prop~qSV-1, data = prop_qsv))
}))

qsv_lm_sacc_df <- cbind(cross_df_sacc,qsv_lm_sacc) 

qsv_lm_sacc_df$p.bonf <- p.adjust(qsv_lm_sacc_df$p.value, 'bonf')
qsv_lm_sacc_df$p.bonf.sig <- qsv_lm_sacc_df$p.bonf < 0.05
qsv_lm_sacc_df$p.fdr<- p.adjust(qsv_lm_sacc_df$p.value, 'fdr')

sum(qsv_lm_sacc_df$p.bonf < 0.05)
# [1] 74
head(qsv_lm_sacc_df)

## top40 
qsv_lm_top40_sacc <- bind_rows(map(cross_sacc, function(x){
  cell = x[[1]]
  qsv = x[[2]]
  prop_qsv <- prop_qsv_sacc %>%
    filter(cell_type == cell & qSV_name == qsv) 
  
  tidy(lm(prop_top40~qSV-1, data = prop_qsv))
}))

qsv_lm_top40_sacc_df <- cbind(cross_df_sacc,qsv_lm_top40_sacc) 

qsv_lm_top40_sacc_df$p.bonf <- p.adjust(qsv_lm_top40_sacc_df$p.value, 'bonf')
qsv_lm_top40_sacc_df$p.bonf.sig <- qsv_lm_top40_sacc_df$p.bonf < 0.05
qsv_lm_top40_sacc_df$p.fdr<- p.adjust(qsv_lm_top40_sacc_df$p.value, 'fdr')

sum(qsv_lm_top40_sacc_df$p.bonf < 0.05)
# [1] 83
head(qsv_lm_top40_sacc_df)

## Amyg
cross_amyg <- cross2(cells_amyg, qsv_names)

cross_df_amyg <- as.data.frame(do.call("rbind", cross_amyg))
colnames(cross_df_amyg) <- c("cell_type","qSV_name")
cross_df_amyg$cell_type <- factor(cross_df_amyg$cell_type, levels = unique(cross_df_amyg$cell_type))
cross_df_amyg$qSV_name <- factor(cross_df_amyg$qSV_name, levels = unique(cross_df_amyg$qSV_name))

## All genes
qsv_lm_amyg <- bind_rows(map(cross_amyg, function(x){
  cell = x[[1]]
  qsv = x[[2]]
  prop_qsv <- prop_qsv_amyg %>%
    filter(cell_type == cell & qSV_name == qsv) 
  
  tidy(lm(prop~qSV-1, data = prop_qsv))
}))

qsv_lm_amyg_df <- cbind(cross_df_amyg, qsv_lm_amyg) 

qsv_lm_amyg_df$p.bonf <- p.adjust(qsv_lm_amyg_df$p.value, 'bonf')
qsv_lm_amyg_df$p.bonf.sig <- qsv_lm_amyg_df$p.bonf < 0.05
qsv_lm_amyg_df$p.fdr<- p.adjust(qsv_lm_amyg_df$p.value, 'fdr')

sum(is.na(qsv_lm_amyg_df$p.value))
# [1] 26
## Replace NA with FALSE 
qsv_lm_amyg_df$p.bonf.sig[is.na(qsv_lm_amyg_df$p.bonf.sig)] <- FALSE
sum(qsv_lm_amyg_df$p.bonf.sig)
# [1] 69

## top40
qsv_lm_top40_amyg <- bind_rows(map(cross_amyg, function(x){
  cell = x[[1]]
  qsv = x[[2]]
  prop_qsv <- prop_qsv_amyg %>%
    filter(cell_type == cell & qSV_name == qsv) 
  
  tidy(lm(prop_top40~qSV-1, data = prop_qsv))
}))

qsv_lm_top40_amyg_df <- cbind(cross_df_amyg, qsv_lm_top40_amyg) 

qsv_lm_top40_amyg_df$p.bonf <- p.adjust(qsv_lm_top40_amyg_df$p.value, 'bonf')
qsv_lm_top40_amyg_df$p.bonf.sig <- qsv_lm_top40_amyg_df$p.bonf < 0.05
qsv_lm_top40_amyg_df$p.fdr<- p.adjust(qsv_lm_top40_amyg_df$p.value, 'fdr')

sum(qsv_lm_top40_amyg_df$p.bonf.sig)
# [1] 85

#### Create scatter plots ####
scatter_qsv_sacc <- prop_qsv_sacc %>% 
  left_join(qsv_lm_sacc_df %>% select(cell_type, qSV_name, p.bonf.sig)) %>% 
  ggplot(aes(qSV, prop, color= p.bonf.sig)) +
  geom_point(size = .4, alpha = .25) +
  facet_grid(cell_type~qSV_name, scales = "free")+
  theme_bw(base_size = 10)+
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'red'))

ggsave(filename = "plots/cellType_qSV_sACC.png", plot = scatter_qsv_sacc, width = 26, height = 10)

scatter_qsv_top40_sacc <- prop_qsv_sacc %>% 
  left_join(qsv_lm_top40_sacc_df %>% select(cell_type, qSV_name, p.bonf.sig)) %>% 
  ggplot(aes(qSV, prop_top40, color= p.bonf.sig)) +
  geom_point(size = .4, alpha = .25) +
  facet_grid(cell_type~qSV_name, scales = "free")+
  theme_bw(base_size = 10)+
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'red'))

ggsave(filename = "plots/cellType_qSV_top40_sACC.png", plot = scatter_qsv_top40_sacc, width = 26, height = 10)

scatter_qsv_amyg <- prop_qsv_amyg %>% 
  left_join(qsv_lm_amyg_df %>% select(cell_type, qSV_name, p.bonf.sig)) %>% 
  ggplot(aes(qSV, prop, color= p.bonf.sig)) +
  geom_point(size = .5, alpha = .25) +
  facet_grid(cell_type~qSV_name, scales = "free")+
  theme_bw(base_size = 10)+
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'red'))

ggsave(filename = "plots/cellType_qSV_amyg.png", plot = scatter_qsv_amyg, width = 26, height = 10)

scatter_qsv_top40_amyg <- prop_qsv_amyg %>% 
  left_join(qsv_lm_top40_amyg_df %>% select(cell_type, qSV_name, p.bonf.sig)) %>% 
  ggplot(aes(qSV, prop_top40, color= p.bonf.sig)) +
  geom_point(size = .5, alpha = .25) +
  facet_grid(cell_type~qSV_name, scales = "free")+
  theme_bw(base_size = 10)+
  scale_color_manual(values = c('FALSE' = 'black', 'TRUE' = 'red'))

ggsave(filename = "plots/cellType_qSV_top40_amyg.png", plot = scatter_qsv_top40_amyg, width = 26, height = 10)
