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
load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)
load(here("differential_expression","data","qSVA_MDD_gene_sex_DEresults.rda"), verbose = TRUE)
outGene_sex = outGene

## Annotate Risk allels ##
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
signifGene <- map_depth(outGene_sex, 3, get_signif, colname = "common_gene_id", return_unique = TRUE)
signifFC <- map_depth(outGene_sex, 3, get_signif, colname = "logFC")
signifSymb <- map_depth(outGene_sex, 3, get_signif, colname = "common_gene_symbol")

## Summarize counts
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
  separate(data, into = c("BrainRegion", "PrimaryDx", "Sex"), sep = "\\.")

model_counts[,1:7]

write_csv(model_counts, here("differential_expression","data","model_counts_gene_sex.csv"))

## Summary table for n signif genes
gene_counts <- model_counts %>%
  dplyr::select(BrainRegion, PrimaryDx, Sex, n_genes) %>%
  pivot_wider(values_from = "n_genes", names_from = "Sex") %>%
  rename(`Control vs.` = PrimaryDx)

write_csv(gene_counts, here("differential_expression","data","model_counts_gene_sex_compare.csv"))

##### Overlap Between Models ####
signifGene_flat <- map(signifGene,my_flatten)

#### Four Way Venns ####
## 1 Venn for each Feature + Model

walk2(signifGene_flat, names(signifGene_flat), function(r_data, r_name){
    filename = paste0("venn_sex_", r_name,'.tiff')
    message(filename)
    
    names(r_data) <- gsub("\\.", " ", names(r_data))
    venn.diagram(r_data, here("differential_expression","plots", filename), 
                 disable.logging = TRUE, 
                 # main = f_name,  
                 # sub= paste("model:", m_name)),
                 fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))
  })

#### Male vs. Female t-stats ####
outGene_sex_compare <- map_depth(outGene_sex, 2, do.call, what = "cbind")
head(outGene_sex_compare$sacc$MDD)

sig_colors <- c(RColorBrewer::brewer.pal(4, "Set1"),"dark grey")
names(sig_colors) <- c("sig_Both", "sig_M", "sig_F", "sig_ALL", "None")

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
  walk2(region_df, names(region_df), 
        ~model_compare_plot(.x, 
                            title = paste(region_name,.y,"vs. Control"),
                            save_png = TRUE,
                            prefix = paste("t-stat_MvF", region_name, .y, sep = "_")))
})
dev.off()

#### Compare Male-only and Female - only with All data
## load all data outGene
load(here("differential_expression","data","qSVA_MDD_gene_DEresults.rda"), verbose = TRUE)
names(outGene$sep)

all(rownames(outGene_sex$amyg$MDD$M) == rownames(outGene$sep$amyg$MDD))

source(here("differential_expression","code","utils.R"))
test <- compare_tstat_plot(DE_out_x = outGene$sep$amyg$MDD,
                           DE_out_y = outGene_sex$amyg$MDD$M, 
                           name_x = "M", 
                           name_y = "All")

ggsave(test, filename = here("differential_expression","plots","test.png"))

pdf(here("differential_expression","plots","sex_vs_All_tstat_scatter.pdf"))

pwalk(list(all_region = outGene$sep, sex_region = outGene_sex, name_region = c("Amygdala", "sACC")), 
      function(all_region, sex_region, name_region){
        pwalk(list(all_dx = all_region, sex_dx = sex_region, name_dx = names(all_region)), 
              function(all_dx, sex_dx, name_dx){
                title1 = paste(name_region, "-", name_dx, "vs. Control,")
                
                walk2(sex_dx, names(sex_dx),
                      ~print(compare_tstat_plot(DE_out_x = all_dx,
                                          DE_out_y = .x, 
                                          name_x = "All", 
                                          name_y = .y,
                                          title = paste(title1, .y))))
              })
      })

dev.off()



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

