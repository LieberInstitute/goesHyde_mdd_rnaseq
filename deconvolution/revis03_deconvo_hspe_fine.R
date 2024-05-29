
library("SummarizedExperiment")
library("SingleCellExperiment")
library("hspe")
library("tidyverse")
library("DeconvoBuddies")
library("here")
library("sessioninfo")

#### Load Data ####
## Load rse_gene data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

rse_gene_region <-  map(c(sacc = "sACC", amy = "Amygdala"), ~rse_gene[,rse_gene$BrainRegion == .x])
map_int(rse_gene_region, ncol)
# sacc  amy 
# 551  540 

## sce Data
load("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/deconvolution_bsp2/data/sce_pan.v2.Rdata", verbose = TRUE)

sce_region <- map(c(sacc = "sacc", amy = "amy"), function(r){
  message(Sys.time(), " - Pseudobulk ", r)
  
  sce <- sce_pan[,sce_pan$region == r]
  sce$cellType <- droplevels(sce$cellType)
  
  sce_pb <- spatialLIBD::registration_pseudobulk(
    sce,
    var_registration = "cellType",
    var_sample_id = "donor",
    covars = NULL,
    min_ncells = 10,
    pseudobulk_rds_file = NULL
  )
  
  return(sce_pb)
})

rm(sce_pan)

map(sce_region, ~table(.x$cellType))
# $sacc
# 
# sacc_Astro_A sacc_Astro_B sacc_Excit_A sacc_Excit_B sacc_Excit_C sacc_Excit_D sacc_Excit_E sacc_Excit_F sacc_Excit_G 
# 4            4            5            5            5            5            5            5            1 
# sacc_Inhib_A sacc_Inhib_B sacc_Inhib_C sacc_Inhib_D sacc_Inhib_E sacc_Inhib_F sacc_Inhib_G sacc_Inhib_H sacc_Inhib_I 
# 5            5            5            5            5            5            4            5            2 
# sacc_Inhib_J sacc_Inhib_K   sacc_Micro sacc_Oligo_A sacc_Oligo_B     sacc_OPC 
# 1            1            4            4            1            4 
# 
# $amy
# 
# amy_Astro_A amy_Astro_B    amy_Endo amy_Excit_A amy_Excit_B amy_Excit_C amy_Inhib_A amy_Inhib_B amy_Inhib_C amy_Inhib_D 
# 5           3           1           4           1           1           2           5           5           5 
# amy_Inhib_E amy_Inhib_F amy_Inhib_G amy_Inhib_H   amy_Micro   amy_Mural   amy_Oligo     amy_OPC   amy_Tcell 
# 1           4           1           1           5           1           5           5           1 

## marker stats
load(here('deconvolution', "data", "revis", "marker_stats_fine.Rdata"), verbose = TRUE)

marker_stats_filter <- map(marker_stats_fine, ~.x |>
                             dplyr::filter(gene %in% rownames(rse_gene) & marker_gene))

#### Run hspe ####
est_prop_hspe <- pmap(list(rse = rse_gene_region, sce = sce_region, markers = marker_stats_filter), 
                        function(rse, sce, markers){
                          
                          pure_samples = rafalib::splitit(sce$cellType)
                          
                          common_genes <- intersect(rowData(sce)$gene_id, rowData(rse)$ensemblID)
                          message("common genes: ", length(common_genes))
                          
                          mixture_samples = t(log2(assays(rse)$counts+1)[common_genes,])
                          reference_samples = t(assays(sce)$logcounts[common_genes,])
                          
                          stopifnot(ncol(mixture_samples) == ncol(reference_samples))
                          
                          marker_genes <- purrr::map(rafalib::splitit(markers$cellType.target), ~markers$gene[.x])
                          marker_genes <- marker_genes[names(pure_samples)]
                          
                          message(Sys.time(), " - run hspe")
                          est_prop = hspe(Y = mixture_samples,
                                               reference = reference_samples,
                                               pure_samples = pure_samples,
                                               markers = marker_genes,
                                               seed = 10524)
                          return(est_prop)
                        })

## Add long data and save
# map(est_prop_hspe, ~round(colMeans(.x$bulk.props), 3))

## save
save(est_prop_hspe, file = here("deconvolution","data" ,"revis","est_prop_hspe_fine.Rdata"))

# slurmjobs::job_single('revis02_deconvo_hspe_fine', create_shell = TRUE, memory = '25G', command = "Rscript revis03_deconvo_hspe_fine.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
