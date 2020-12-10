
library(SummarizedExperiment)
library(SingleCellExperiment)
library(jaffelab)
library(here)
library(MuSiC)
library(Biobase)
library(xbioc)
library(purrr)
library(tidyverse)
library(reshape2)

#### Load Data ####
## Load rse_gene data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sce Data
load(here("deconvolution","data","sce.sacc_filtered.Rdata"), verbose = TRUE)
load(here("deconvolution","data","sce.amyg_filtered.Rdata"), verbose = TRUE)

## marker data
load(here("deconvolution","data","marker_stats.Rdata"), verbose = TRUE)
top_n <- 5

marker_genes <- map(marker_stats, ~.x %>%
                      arrange(-ratio) %>%
                      filter(rank_ratio <= top_n) %>%
                      arrange(cellType.target, rank_ratio) %>% 
                      pull(gene))

map_int(marker_genes, length)
map_int(marker_genes, ~sum(.x %in% rownames(rse_gene)))


exp <- list(sacc = list(bulk = rse_gene[,rse_gene$BrainRegion == "sACC"],
                        sce = sce.sacc), 
            amyg = list(bulk = rse_gene[,rse_gene$BrainRegion == "Amygdala"],
                        sce = sce.amyg)
)
map(exp, ~map_int(.x, ncol))
# $sacc
# bulk  sce 
# 551 7004 
# 
# $amyg
# bulk  sce 
# 540 6582 
exp_set <- map(exp, ~list(bulk = ExpressionSet(assayData = assays(.x$bulk)$counts),
                          sce = ExpressionSet(assayData = as.matrix(assays(.x$sce)$counts),
                                              phenoData=AnnotatedDataFrame(
                                                as.data.frame(colData(.x$sce))[c("cellType","cellType.Broad", "uniqueID","donor")])
                                              )
                          )
               )
      

#### estimate cell type props ####
ct <- list(broad = "cellType.Broad", specific = "cellType")

est_prop <- map2(marker_genes, names(marker_genes), function(mg, n){
    est_prop <- music_prop(bulk.eset = exp_set[[ss(n, "_")]]$bulk,
             sc.eset = exp_set[[ss(n, "_")]]$sce,
             clusters = ct[[ss(n, "_",2)]],
             samples = 'donor',
             markers = mg)
    return(est_prop)
} 
)

map(est_prop, ~round(colMeans(.x$Est.prop.weighted),3))
# $sacc_broad
# Oligo Micro Astro Inhib   OPC Excit 
# 0.156 0.052 0.552 0.211 0.000 0.028 
# 
# $amyg_broad
# Inhib Oligo Astro Excit Micro   OPC 
# 0.192 0.169 0.474 0.085 0.072 0.008 
# 
# $sacc_specific
# Oligo   Micro   Astro Inhib.2 Inhib.1     OPC Excit.3 Excit.1 Excit.2 Excit.4 
# 0.220   0.080   0.595   0.011   0.035   0.003   0.017   0.038   0.001   0.000 
# 
# $amyg_specific
# Inhib.2   Oligo   Astro Excit.1   Micro     OPC Inhib.3 Inhib.5 Inhib.1 
# 0.023   0.296   0.348   0.051   0.189   0.003   0.001   0.053   0.035

save(est_prop, file= here("deconvolution","data",paste0("est_prop_top", top_n,".Rdata")))

pd <- as.data.frame(colData(rse_gene))

pd2 <- pd %>%
  select(PrimaryDx, ERCCsumLogErr, grep("snpPC",colnames(pd))) %>%
  rownames_to_column("sample")

est_prop_long <- map(est_prop, ~melt(.x$Est.prop.weighted) %>%
                       rename(sample = Var1, cell_type = Var2, prop = value) %>% 
                       arrange(cell_type) %>%
                       left_join(pd2, by = "sample"))

save(est_prop_long, file = here("deconvolution","data",paste0("est_prop_top", top_n,"_long.Rdata")))

# sgejobs::job_single('music_deconvo', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript music_deconvo.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


