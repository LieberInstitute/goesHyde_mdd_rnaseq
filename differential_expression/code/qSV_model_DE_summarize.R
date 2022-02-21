library(jaffelab)
library(SummarizedExperiment)
library(tidyverse)
library(sessioninfo)
library(here)
library(EnhancedVolcano)
library(UpSetR)
library(VariantAnnotation)
library(pheatmap)
library(VennDiagram)

source(here("differential_expression","code", "utils.R"))

#### Load Data ####
data_type <- c("gene","exon","jxn","tx")
names(data_type) <- data_type

paths <- map(data_type, ~here("differential_expression","data",paste0("qSVA_MDD_", .x, "_DEresults.rda")))

allOut <- lapply(paths, function(x) get(load(x, verbose = TRUE)))

load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)

## Annotate Risk allels ##
risk_vcf <- readVcf(here("eqtl", "data", "risk_snps", "LIBD_maf01_gwas_MDD.vcf.gz"))
mdd_snps <- read.delim(here("eqtl", "data", "risk_snps", "PGC_MDD_genome-wide_significant_Jan2022.txt")) %>%
  filter(!is.na(bp))
## find genes 1KB from risk snp
risk_gr <- GRanges(seqnames = mdd_snps$chr, IRanges(mdd_snps$bp, width = 1), feature_id = mdd_snps$markername)
oo <- findOverlaps(risk_gr, rowRanges(rse_gene), maxgap = 100000)
risk_genes <- unique(rowRanges(rse_gene)[subjectHits(oo),])
length(risk_genes$gencodeID)
# [1] 538

#### Summarize Counts ####
## Extract lists
signifFeat <- map_depth(allOut, 4, get_signif)
signifGene <- map_depth(allOut, 4, get_signif, colname = "common_gene_id", return_unique = TRUE)
signifGeneRep <- map_depth(allOut, 4, get_signif, colname = "common_gene_id", return_unique = FALSE) ## non-unique genes to tabulate
signifFC <- map_depth(allOut, 4, get_signif, colname = "logFC")
signifSymb <- map_depth(allOut, 4, get_signif, colname = "common_gene_symbol")

## Summarize counts
signifFeat_n <- map_depth(signifFeat, 4, length) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(n_features = V1)

signifGene_n <- map_depth(signifGene, 4, length) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(n_genes = V1)

signifRiskGene_n <- map_depth(signifGene, 4, ~sum(.x %in% risk_genes$gencodeID)) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(n_risk_genes = V1)

signifGene_n_up <- map_depth(signifFC, 4, ~sum(.x > 0)) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(Up = V1)

signifGene_n_down <- map_depth(signifFC, 4, ~sum(.x < 0)) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(Down = V1)

signifSymbol <- map_depth(signifSymb, 4, ~paste(unique(.x)[order(unique(.x))], collapse = ", ")) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(Genes = V1)

signifRiskSymbol <- map_depth(signifSymb, 4, ~paste(unique(.x[.x %in% risk_genes$Symbol]), collapse = ", ")) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rename(RiskGenes = V1)

model_counts <- do.call("cbind",list(signifFeat_n, 
                                     signifGene_n, 
                                     signifGene_n_up, 
                                     signifGene_n_down,
                                     signifRiskGene_n, 
                                     signifRiskSymbol, 
                                     signifSymbol)) %>%
  rownames_to_column("data") %>%
  separate(data, into = c("Feature", "Model", "Region", "Dx"), sep = "\\.")

model_counts[,1:9]

write_csv(model_counts, here("differential_expression","data","summary_tables","model_counts.csv"))
# model_counts <- read_csv(here("differential_expression","data","summary_tables","model_counts.csv"))

## Get Number of genes across features
gene_counts <- model_counts %>% filter(Model == "sep") %>%
  dplyr::select(Region, Dx, Feature, n_genes) %>%
  arrange(Region, Dx)

write.csv(gene_counts, file =  here("differential_expression","data","summary_tables","sep_gene_counts.csv"))

unique_gene_counts <- model_counts %>% filter(Model == "sep") %>%
  group_by(Region, Dx) %>%
  summarize(non_unique_genes = sum(n_genes)) %>%
  arrange(Region, Dx)


## Unique genes by 
# signifGene$gene$sep$amyg$MDD
signifGene_flip <-map(transpose(transpose(signifGene)$sep), transpose)

unique_genes <- my_flatten(map_depth(signifGene_flip, 2, ~unique(unlist(.x))))

(unique_gene_counts <- unique_gene_counts %>%
  left_join(t(map_dfc(unique_genes, length)) %>%
  as.data.frame() %>%
  rownames_to_column("group") %>%
  separate(group, into = c("Region", "Dx")) %>%
  rename(unique_genes = V1)))

## Updates in GO_enrichment
write.csv(unique_gene_counts, file =  here("differential_expression","data","summary_tables","sep_gene_counts_unique.csv"))

## Venn Diagrams
walk2(signifGene_flip, names(signifGene_flip), function(region_data, region_name){
  walk2(region_data, names(region_data), function(dx_data, dx_name){
    
    filename = paste0("venn_unique_genes_", region_name, "_", dx_name,'.tiff')
    message(filename)
    
    venn.diagram(dx_data, here("differential_expression","plots", filename), 
                 disable.logging = TRUE, 
                 main = "Overlapping Genes Between Features",
                 sub = paste("Region:", region_name, "Dx:", dx_name),
                 fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))
  })
})


venn.diagram(unique_genes, here("differential_expression","plots", "venn_unique_genes_ALL.png"), 
             disable.logging = TRUE, 
             main = "Unique Genes",
             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))

#### heatmaps ####
summary_pheat <- function(count){
  
  heat_data <- model_counts %>%
    mutate(name = paste(Model, Region, Dx)) %>%
    dplyr::select(all_of(c("name", "Feature", count))) %>%
    pivot_wider(names_from = "Feature", values_from = count) %>%
    column_to_rownames("name") 
  
  png(here("differential_expression","plots",paste0("heat_",count,".png")), height = 1000, width = 800)
  pheatmap(heat_data,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           display_numbers = TRUE,
           number_format = "%i",
           fontsize = 20,
           gaps_row = c(2, 4, 6),
           main = count)
  dev.off()
}

map(c("n_genes", "n_features", "Up", "Down", "n_risk_genes"), summary_pheat)

#### How may significant features are from one genes? ###
signifGeneRep <- map_depth(signifGeneRep, 4, table)
map_depth(signifGeneRep, 4, max)
map_depth(signifGeneRep, 4, ~sum(.x > 1)/length(.x))

signifGeneRep2 <- transpose(signifGeneRep)$sep
signifGeneRep_flat <- my_flatten(signifGeneRep2)
names(signifGeneRep_flat)

## Do heat map?
signifGeneRep_df <- map2(signifGeneRep_flat, names(signifGeneRep_flat),function(t, name){
  df <- as.data.frame(t)
  df$Set <- name
  return(df)
})

signifGeneRep_long <- do.call("rbind", signifGeneRep_df)
rownames(signifGeneRep_long) <- NULL
signifGeneRep_long %>% arrange(-Freq) %>% head

signifGeneRep_tab <- signifGeneRep_long %>%
  filter(Freq > 1) %>%
  pivot_wider(values_from = "Freq", names_from = "Set") %>%
  column_to_rownames("Var1")

corner(signifGeneRep_tab)
dim(signifGeneRep_tab)
# [1] 2784   16

signifGeneRep_tab0 <- signifGeneRep_tab
signifGeneRep_tab0[is.na(signifGeneRep_tab0)]  <- 0
# signifGeneRep_tab0[signifGeneRep_tab0 > 10 ]  <- 10

png(here("differential_expression","plots","heat_gene_reps.png"), height = 800, width = 550)
pheatmap(signifGeneRep_tab, show_rownames = FALSE, cluster_rows=FALSE, cluster_cols=FALSE)
# pheatmap(signifGeneRep_tab0, show_rownames = FALSE)
dev.off()

## Density?
gene_rep_denisty <- signifGeneRep_long %>%
  filter(!grepl("gene", Set), Freq > 1) %>%
  ggplot(aes(x = Freq, color = Set)) +
  geom_density()

ggsave(gene_rep_denisty, filename = here("differential_expression","plots","gene_rep_density.png"))

signifGeneRep_long %>%
  filter(!grepl("gene", Set), Freq > 1) %>%
  count(Set)

##### Overlap Between Models ####
signifFeat_flat <- map(signifFeat,my_flatten)

#### Using map causes bug?
pdf(here("differential_expression","plots",paste0("upset_","gene",".pdf")))
upset(fromList(signifFeat_flat$gene), order.by = "freq", nset = 8)
dev.off()

pdf(here("differential_expression","plots",paste0("upset_exon.pdf")))
upset(fromList(signifFeat_flat$exon), order.by = "freq", nset = 8)
dev.off()

pdf(here("differential_expression","plots",paste0("upset_jxn.pdf")))
upset(fromList(signifFeat_flat$jxn), order.by = "freq", nset = 8)
dev.off()

pdf(here("differential_expression","plots",paste0("upset_tx.pdf")))
upset(fromList(signifFeat_flat$tx), order.by = "freq", nset = 8)
dev.off()


##only sep
walk2(signifFeat, names(signifFeat), function(sf, name){
  pdf(here("differential_expression","plots",paste0("upset_",name,"_sep",".pdf")))
  print(upset(fromList(my_flatten(sf$sep)), order.by = "freq", nset = 4))
  dev.off()
})

signifGene_flat <- my_flatten(signifGene)
pdf(here("differential_expression","plots",paste0("upset_ALL.pdf")))
upset(fromList(signifGene_flat), order.by = "freq", nset = sum(model_counts$n_genes != 0))
dev.off()

#### Model vs. Model ####
## Compare Model Intersects
# names(signifFeat$gene$sep$amyg$MDD)
signifFeat_model <- map(signifFeat, function(sg) map(transpose(sg), transpose))
signifFeat_region <- map(signifFeat, function(sg) map(sg, transpose))
# names(signifFeat_region$gene$sep$MDD)

get_venn_values <- function(my_sets){
  my_values <- c(length(setdiff(my_sets[[1]], my_sets[[2]])),
                    length(intersect(my_sets[[1]], my_sets[[2]])),
                    length(setdiff(my_sets[[2]], my_sets[[1]])))
  names(my_values) <- c(names(my_sets)[[1]],
                      "Intersect", names(my_sets)[[2]])
  return(my_values)
}

model_venn_values <- as.data.frame(map(signifFeat_model, function(sg) map(sg, ~map(.x, get_venn_values)))) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("cat") %>%
  pivot_longer(!cat)%>%
  separate(cat, into = c("Feature", "Region", "Dx"), sep = "\\.")

model_venn_values$Feature <- factor(model_venn_values$Feature, levels = c("gene", "exon", "jxn", "tx"))
model_venn_values$name <- factor(model_venn_values$name,  levels = c("sep", "Intersect", "sep_cf"))

model_venn_tiles <- ggplot(model_venn_values, aes(x = name, y = Feature, fill = value)) +
  geom_tile(color = "lightgrey") +
  geom_text(aes(label = value), color = "white") +
  labs(x = NULL, y = NULL)+
  facet_grid(Region~Dx)
  
ggsave(model_venn_tiles, filename = here("differential_expression","plots", "venn_tiles_models.png"))

## Region venn values
region_venn_values <- as.data.frame(map(signifFeat_region, function(sg) map(sg, ~map(.x, get_venn_values)))) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("cat") %>%
  pivot_longer(!cat) %>%
  separate(cat, into = c("Feature", "Model", "Dx"), sep = "\\.") 
# %>%
#   dplyr::filter(grepl("sep", Model))

region_venn_values$Feature <- factor(region_venn_values$Feature, levels = c("gene", "exon", "jxn", "tx"))

region_venn_tiles <- ggplot(region_venn_values, aes(x = name, y = Feature, fill = value)) +
  geom_tile(color = "lightgrey") +
  geom_text(aes(label = value), color = "white") +
  labs(x = NULL, y = NULL) +
  facet_grid(Model~Dx)

ggsave(region_venn_tiles, filename = here("differential_expression","plots", "venn_tiles_region.png"))

#### Four Way Venns ####
## 1 Venn for each Feature + Model
signifFeat_venn_data <- map_depth(signifFeat, 2, my_flatten)
map_depth(signifFeat_venn_data, 2, names)

walk2(signifFeat_venn_data, names(signifFeat_venn_data), function(f_data, f_name){
  walk2(f_data, names(f_data), function(m_data, m_name){
    filename = paste0("venn_",f_name, "_", m_name,'.jpg')
    message(filename)
    
    names(m_data) <- gsub("\\.", " ", names(m_data))
    venn.diagram(m_data, here("differential_expression","plots", filename), 
                 disable.logging = TRUE, 
                 # main = f_name,  
                 # sub= paste("model:", m_name)),
                 fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))
  })
})
  
## Graph t-stats
allOut_model <- map(allOut, function(x) map(transpose(x), transpose))
names(allOut_model$gene$amyg$MDD)

head(do.call("cbind", allOut_model$gene$amyg$MDD))

allOut_model_compare <- map_depth(allOut_model, 3, do.call, what = "cbind")

sig_colors <- c(RColorBrewer::brewer.pal(3, "Set1"),"dark grey")
names(sig_colors) <- c("sig_Both", "sig_no-deconvo", "sig_deconvo", "None")

model_compare_plot <- function(model_df, title){

  t_cor <- model_df %>%
    summarize(cor = cor(sep.t, sep_cf.t, method = "spearman")) %>%
    mutate(cor_anno = paste0("rho==", format(round(cor, 2), nsmall = 2)))
  
  model_df <- model_df %>% 
    mutate(Signif = case_when(sep.adj.P.Val < 0.05 & sep_cf.adj.P.Val < 0.05 ~"sig_Both",
                              sep.adj.P.Val < 0.05 ~ "sig_no-deconvo",
                              sep_cf.adj.P.Val < 0.05 ~ "sig_deconvo",
                              TRUE ~ "None"))
  
  model_plot <- ggplot(model_df, aes(x = sep.t, y = sep_cf.t, color = Signif)) +
    geom_point(alpha = 0.5, size = 0.5) +
    labs(x = "t-stats DE Model", y = "t-stats DE Model with Cell Fractions", 
         title = title, subtitle = t_cor$cor_anno) +
    scale_color_manual(values = sig_colors)+
    theme_bw()
  
  message(title)  
  print(model_plot)
}

pdf(here("differential_expression","plots","qSV_model_scatter.pdf"))
walk2(allOut_model_compare, names(allOut_model_compare), function(feat_df, feat_name){
  walk2(feat_df, names(feat_df), function(region_df, region_name){
    walk2(region_df, names(region_df), ~model_compare_plot(.x, title = paste(feat_name, region_name,.y,"vs. Control")))
  })
})
dev.off()

## Compare Dx or Region?

#### Volcano Plots ####
## need to fix Tx output!
pdf(here("differential_expression","plots","qSV_model_volcano.pdf"))
map2(allOut, c("Gene", "Exon", "Jxn", "Tx"), function(featOut, featName){
  map2(featOut, names(featOut), function(modelOut, modelName){
    map2(modelOut, names(modelOut), function(regionOut, regionName){
      map2(regionOut, names(regionOut),
           ~EnhancedVolcano(.x,
                            pCutoff = 0.05,
                            y = "adj.P.Val",
                            x = "logFC",
                            lab = .x$Symbol,
                            title = paste(featName,"-", regionName),
                            subtitle = paste(modelName, "model,", .y,"vs. Control"))
      )
    })
  })
})

dev.off()

#### Boxplots ####

# ranks_to_check <- c(1:50, 201:250, 401:450, 601:650, 801:850, 1001:1050)
ranks_to_check <- c(1:50)
length(ranks_to_check)

outGene_test <- outGene$sep$amyg$MDD %>%
  arrange(adj.P.Val) %>%
  mutate(row_num = row_number())


#sgejobs::job_single('qSV_model_DE_summarize', queue = 'bluejay', create_shell = TRUE, memory = '80G', command = "Rscript qSV_model_DE_summarize.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

