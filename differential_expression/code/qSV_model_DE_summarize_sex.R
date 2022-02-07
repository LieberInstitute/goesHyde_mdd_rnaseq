library(jaffelab)
library(SummarizedExperiment)
library(tidyverse)
library(sessioninfo)
library(here)
# library(EnhancedVolcano)
# library(UpSetR)
library(VariantAnnotation)
# library(pheatmap)
library(VennDiagram)

source(here("differential_expression","code","utils.R"))

#### Load Data ####
load(here("differential_expression","data","qSVA_MDD_gene_sex_DEresults.rda"), verbose = TRUE)

# ## Annotate Risk allels ##
risk_vcf <- readVcf(here("eqtl", "data", "risk_snps", "LIBD_maf01_gwas_MDD.vcf.gz"))
mdd_snps <- read.delim(here("eqtl", "data", "risk_snps", "PGC_MDD_genome-wide_significant_Jan2022.txt")) %>%
  filter(!is.na(bp))

## find genes 1KB from risk snp
risk_gr <- GRanges(seqnames = paste0(mdd_snps$chr), IRanges(mdd_snps$bp, width = 1), feature_id = mdd_snps$markername)
oo <- findOverlaps(risk_gr, rowRanges(rse_gene), maxgap = 100000)
risk_genes <- unique(rowRanges(rse_gene)[subjectHits(oo),])
length(risk_genes$gencodeID)
# [1] 538

#### Summarize Counts ####
## Extract lists
signifFeat <- map_depth(outGene, 3, get_signif)
signifGene <- map_depth(outGene, 3, get_signif, colname = "common_gene_id", return_unique = TRUE)
signifFC <- map_depth(outGene, 3, get_signif, colname = "logFC")
signifSymb <- map_depth(outGene, 3, get_signif, colname = "common_gene_symbol")

## Summarize counts
# signifFeat_n <- map_depth(signifFeat, 3, length) %>%
#   as.data.frame() %>%
#   t() %>%
#   as.data.frame() %>%
#   rename(n_features = V1)

signifGene_n <- map_depth(signifGene, 3, length) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(n_genes = V1)

signifRiskGene_n <- map_depth(signifGene, 3, ~sum(.x %in% risk_genes$gencodeID)) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(n_risk_genes = V1)

signifGene_n_up <- map_depth(signifFC, 3, ~sum(.x > 0)) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(Up = V1)

signifGene_n_down <- map_depth(signifFC, 3, ~sum(.x < 0)) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(Down = V1)

signifSymbol <- map_depth(signifSymb, 3, ~paste(unique(.x)[order(unique(.x))], collapse = ", ")) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(Genes = V1)

signifRiskSymbol <- map_depth(signifSymb, 3, ~paste(unique(.x[.x %in% risk_genes$Symbol]), collapse = ", ")) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(RiskGenes = V1)

model_counts <- do.call("cbind",list(signifGene_n, signifGene_n_up, signifGene_n_down,signifRiskGene_n, signifRiskSymbol, signifSymbol)) %>%
  rownames_to_column("data") %>%
  separate(data, into = c("Region", "Dx", "Sex"), sep = "\\.")

model_counts[,1:7]

write_csv(model_counts, here("differential_expression","data","model_counts_gene_sex.csv"))

##### Overlap Between Models ####
my_flatten <- function (x, use.names = TRUE, classes = "ANY") 
{
  #' Source taken from rlist::list.flatten
  len <- sum(rapply(x, function(x) 1L, classes = classes))
  y <- vector("list", len)
  i <- 0L
  items <- rapply(x, function(x) {
    i <<- i + 1L
    y[[i]] <<- x
    TRUE
  }, classes = classes)
  if (use.names && !is.null(nm <- names(items))) 
    names(y) <- nm
  y
}

signifGene_flat <- map(signifGene,my_flatten)


#### Four Way Venns ####
## 1 Venn for each Feature + Model

walk2(signifGene_flat, names(signifGene_flat), function(r_data, r_name){
    filename = paste0("venn_sex_", r_name,'.jpg')
    message(filename)
    
    names(r_data) <- gsub("\\.", " ", names(r_data))
    venn.diagram(r_data, here("differential_expression","plots", filename), 
                 disable.logging = TRUE, 
                 # main = f_name,  
                 # sub= paste("model:", m_name)),
                 fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))
  })

#### Male vs. Female t-stats ####
outGene_sex_compare <- map_depth(outGene, 2, do.call, what = "cbind")
head(outGene_sex_compare$sacc$MDD)

sig_colors <- c(RColorBrewer::brewer.pal(3, "Set1"),"dark grey")
names(sig_colors) <- c("sig_Both", "sig_M", "sig_F", "None")

model_compare_plot <- function(model_df, title, save_png = FALSE, prefix = ""){

  t_cor <- model_df %>%
    summarize(cor = cor(M.t, F.t, method = "spearman")) %>%
    mutate(cor_anno = paste0("rho = ", format(round(cor, 2), nsmall = 2)))
  
  model_df <- model_df %>% 
    mutate(Signif = case_when(M.adj.P.Val < 0.05 & F.adj.P.Val < 0.05 ~"sig_Both",
                              M.adj.P.Val < 0.05 ~ "sig_M",
                              F.adj.P.Val < 0.05 ~ "sig_F",
                              TRUE ~ "None"))
  
  model_plot <- ggplot(model_df, aes(x = M.t, y = F.t, color = Signif)) +
    geom_point(alpha = 0.7, size = 0.5) +
    labs(x = "t-stats Male", y = "t-stats Female", 
         title = title, subtitle = t_cor$cor_anno, parse = T) +
    scale_color_manual(values = sig_colors)+
    theme_bw()
  
  message(title)  
  
  if(save_png){
    fn = here("differential_expression","plots",paste0(prefix,".png"))
    ggsave(model_plot, filename = fn)              
  }
  print(model_plot)
}

pdf(here("differential_expression","plots","sex_tstat_scatter.pdf"))
  walk2(outGene_sex_compare, c("Amygdala", "sACC"), function(region_df, region_name){
    walk2(region_df, names(region_df), ~model_compare_plot(.x, 
                                                           title = paste(region_name,.y,"vs. Control"),
                                                           save_png = TRUE,
                                                           prefix = paste("t-stat_MvF", region_name, .y, sep = "_")))
  })
dev.off()

#### Compare Male-only with All data
outGene_sex <- outGene
names(outGene_sex$sacc$MDD)

outGene_M <- transpose(map(outGene_sex, transpose))[["M"]]
names(outGene_M)

## load all data outGene
load(here("differential_expression","data","qSVA_MDD_gene_DEresults.rda"), verbose = TRUE)
names(outGene$sep)

all(rownames(outGene_M$amyg$MDD) == rownames(outGene$sep$amyg$MDD))

compare_sex_all <- map2(outGene_M, outGene$sep, function(region_M, region_all){
  map2(region_M, region_all, function(dx_M, dx_all){
    colnames(dx_M) <- paste0("M.", colnames(dx_M))
    return(cbind(dx_all, dx_M))
  })
})
head(compare_sex_all$sacc$MDD)

sig_colors <- c(RColorBrewer::brewer.pal(3, "Set1"),"dark grey")
names(sig_colors) <- c("sig_Both", "sig_M", "sig_all", "None")

M_compare_plot <- function(model_df, title, save_png = FALSE, prefix = ""){
  
  t_cor <- model_df %>%
    summarize(cor = cor(t, M.t, method = "spearman")) %>%
    mutate(cor_anno = paste0("rho = ", format(round(cor, 2), nsmall = 2)))
  
  model_df <- model_df %>% 
    mutate(Signif = case_when(M.adj.P.Val < 0.05 & adj.P.Val < 0.05 ~"sig_Both",
                              M.adj.P.Val < 0.05 ~ "sig_M",
                              adj.P.Val < 0.05 ~ "sig_all",
                              TRUE ~ "None"))
  
  model_plot <- ggplot(model_df, aes(x = t, y = M.t, color = Signif)) +
    geom_point(alpha = 0.7, size = 0.5) +
    labs(x = "t-stats All", y = "t-stats Male-only", 
         title = title, subtitle = t_cor$cor_anno) +
    scale_color_manual(values = sig_colors)+
    theme_bw()
  
  if(save_png){
    fn = here("differential_expression","plots",paste0(prefix,".png"))
    ggsave(model_plot, filename = fn)              
  }
  
  message(title)  
  print(model_plot)
}

pdf(here("differential_expression","plots","sexM_vs_All_tstat_scatter.pdf"))
walk2(compare_sex_all, c("Amygdala", "sACC"), function(region_df, region_name){
  walk2(region_df, names(region_df), ~M_compare_plot(.x, 
                                                     title = paste(region_name,.y,"vs. Control"),
                                                     save_png = TRUE,
                                                     prefix = paste("t-stat_MvAll", region_name, .y, sep = "_")))
})
dev.off()

#### T-stat distribution ####
## compare M and F
load(here("differential_expression","data","qSVA_MDD_gene_M_downsample_DEresults.rda"), verbose = TRUE)
length(outGene_downsample$amyg$MDD)

plot_t_density <- function(DE_main, DE_downsample, title){
  
  DE_main <- map2(DE_main,names(DE_main), ~.x %>% mutate(Sex = .y))
  DE_main <- do.call("rbind", DE_main) 
  
  DE_downsample <- map2(DE_downsample,1:length(DE_downsample), ~.x %>% mutate(rep = .y))
  DE_downsample <- do.call("rbind", DE_downsample) 

  t_density <- ggplot() +
    geom_density(data = DE_downsample, aes(x = t, group = rep), size=0.2, colour=alpha("dark grey", 0.3))+
    geom_density(data = DE_main,  aes(x = t, color = Sex))+
    labs(title = title) +
    theme_bw() +
    theme(text = element_text(size=15))

  return(t_density)
}


pwalk(list(main_region = outGene_sex, ds_region = outGene_downsample, region_name = names(outGene_sex)),
     function(main_region, ds_region, region_name){
       pwalk(list(main_dx = main_region, ds_dx = ds_region, dx_name = names(main_region)),
            function(main_dx, ds_dx, dx_name){
              t_density <- plot_t_density(DE_main = main_dx, DE_downsample = ds_dx, title = paste(region_name, ":", dx_name, "vs. Control"))
              ggsave(t_density, filename = here("differential_expression","plots", paste0("t_stats_distribution_sex-",region_name, "_",dx_name, ".png")))
            })
     })


#sgejobs::job_single('qSV_model_DE_summarize', queue = 'bluejay', create_shell = TRUE, memory = '80G', command = "Rscript qSV_model_DE_summarize.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

