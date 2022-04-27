#####################################
# from /dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/R_scripts/DEG_scripts/qSV_model_analysis_DLPFC_PTSD.R
#####################################

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)
library(edgeR)
library(tidyverse)
library(here)
library(RColorBrewer)

#### load objects ####

#### expression ####
load(here('exprs_cutoff','rse_gene.Rdata'), verbose=TRUE)

rse_gene$Dx = as.factor(rse_gene$PrimaryDx)
rse_gene$Dx = factor(rse_gene$Dx, levels(rse_gene$Dx)[c(2,3,1)] )
levels(rse_gene$Dx)

rd <- as.data.frame(rowData(rse_gene))
pd <- as.data.frame(colData(rse_gene))

## Plots needed
april_plots <- readxl::read_excel(here("eqtl","data","plot_data","Genes_for_boxplot_042622.xlsx"), sheet = "Diff_exp")
april_plots <- april_plots %>% rename(Symbol = symbol) %>% left_join(rd %>% select(Symbol, gencodeID))

## load DE results
load(here("differential_expression","data","qSVA_MDD_gene_DEresults.rda"), verbose=TRUE)
## only need sep
outGene <- outGene$sep
names(outGene)

outGene_april <- map_depth(outGene, 2, ~.x[april_plots$gencodeID,])
map_depth(outGene_april, 2, dim)

## match by Symbol
all(april_plots$symbol %in% rowData(rse_gene)$Symbol)

## load degradation data
load(here("differential_expression","data","degradation_rse_MDDseq_BiPSeq_BothRegions.Rdata"), verbose = TRUE)
## load models
load(here("differential_expression","data","differental_models.Rdata"), verbose = TRUE)
#load qSVs
load(here("differential_expression","data","qSV_mat.Rdata"), verbose = TRUE)


## Get resid expression ##
gRpkm = recount::getRPKM(rse_gene,"Length")
gExprs = as.matrix(log2(gRpkm+1))

gExprs <- map(splitit(rse_gene$BrainRegion), ~gExprs[,.x])
# gExprs <- map2(splitit(rse_gene$BrainRegion), modSep, ~cleaningY(gExprs[,.x], .y, P=3)) ## test non MDD and not BPD

cleaningY(gExprs$Amygdala, modSep$amyg$MDD, P= 3)

gExprs_resid <- map2(gExprs, modSep, function(gExprs_region, models){
  map(models, function(m){
    gExprs_dx <- gExprs_region[,rownames(m)]
    message("samples match:", all(rownames(m) == colnames(gExprs_dx)))
    cleaningY(gExprs_dx, m, P=3)
  })
})
map_depth(gExprs_resid, 2, dim)

pd_split <- map_depth(modSep, 2, ~pd[rownames(.x),])
map_depth(pd_split, 2, corner)

# topGenes_residExp_beeswarm(outGene_april$amyg$MDD, gExprs_resid$amyg$MDD, pd_split$amyg$MDD, "plots/april_test.pdf")
# 

## set up plotting
plot_dir <- here("differential_expression","plots","april_plots")
load(here("data", "MDD_colors.Rdata"), verbose = TRUE)

de_plot_gene <- function(e = gExprs_resid$Amygdala$MDD, gene = "ENSG00000038219.12", pd = pd_split$amyg$MDD, title, subtitle){
  
  pd$PrimaryDx <- droplevels(pd$PrimaryDx)
  dx_levels <- levels(pd$PrimaryDx)
  
  e_subset <- as.data.frame(e[gene,]) %>%
    rownames_to_column("Sample") %>%
    add_column(PrimaryDx = pd$PrimaryDx) %>%
    rename(resid_expression = `e[gene, ]`)
  
  DE_box <- ggplot(data = e_subset,
                   aes(x = PrimaryDx, y = resid_expression, fill = PrimaryDx)) +
    geom_boxplot() +
    scale_fill_manual(values = mdd_Dx_colors[dx_levels]) +
    theme_bw() +
    theme(text = element_text(size = 15)) +
    labs(x = "Primary Dx", y = "Residual Expression", title = title, subtitle = subtitle)
  
  return(DE_box)
}


de_plot <- function(e = gExprs_resid$Amygdala$MDD, og = outGene_april$amyg$MDD, pd = pd_split$amyg$MDD, file){
 og_anno <- og %>% 
   mutate(anno = paste0("FDR = ", round(adj.P.Val,3), ", logFC = ", round(logFC,3), "\nensemblID:",ensemblID)) %>% 
   select(Symbol, gencodeID, anno)
 
 pdf(here(plot_dir,paste0(file,".pdf")))
 pmap(og_anno, function(Symbol, gencodeID, anno){
   # message(Symbol," gene = ",gencodeID)
   print(de_plot_gene(e = e, gene = gencodeID, pd = pd, title = Symbol, subtitle = anno))
 })
 dev.off()
}



de_plot(e = gExprs_resid$Amygdala$MDD, og = outGene_april$amyg$MDD, pd = pd_split$amyg$MDD, file = "AprilPlots_DE_amyg_MDD")
de_plot(e = gExprs_resid$sACC$MDD, og = outGene_april$sacc$MDD, pd = pd_split$sacc$MDD, file = "AprilPlots_DE_sacc_MDD")
de_plot(e = gExprs_resid$Amygdala$BPD, og = outGene_april$amyg$BPD, pd = pd_split$amyg$BPD, file = "AprilPlots_DE_amyg_BPD")
de_plot(e = gExprs_resid$sACC$BPD, og = outGene_april$sacc$BPD, pd = pd_split$sacc$BPD, file = "AprilPlots_DE_sacc_BPD")



ggsave(test, file = here(plot_dir,"test_box.png"))

#sgejobs::job_single('qSV_model_DE_plot', create_shell = TRUE, memory = '80G', command = "Rscript qSV_model_DE_plot.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
