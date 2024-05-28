
library("DeconvoBuddies")
library("SingleCellExperiment")
library("purrr")
library("dplyr")
library("here")
library("sessioninfo")

## sce Data
load("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/deconvolution_bsp2/data/sce_pan.v2.Rdata", verbose = TRUE)

sce_region <- map(c(sacc = "sacc", amy = "amy"), function(r){
  sce <- sce_pan[,sce_pan$region == r]
  sce$cellType <- droplevels(sce$cellType)
  return(sce)
})
rm(sce_pan)

map_int(sce_region, ncol)
# sacc   amy 
# 15323 14006 
map(sce_region, ~table(.x$cellType))
map_int(sce_region, ~length(levels((.x$cellType))))
# sacc  amy 
# 24   19
marker_stats_fine <- map(sce_region, ~get_mean_ratio(.x, 
                               cellType_col = "cellType",
                               gene_ensembl = "gene_id",
                               gene_name = "gene_name") |>
                      mutate(marker_gene = MeanRatio.rank <=25 & 
                               MeanRatio > 1)
)


map(marker_stats_fine, ~.x |>
      dplyr::filter(marker_gene) |>
      dplyr::count(cellType.target) |>
      arrange(n))

# $sacc
# cellType.target     n
#   1 sacc_Oligo_B       21
# 
# $amy
# cellType.target     n
#   1 amy_Astro_B         5

save(marker_stats_fine, file = here('deconvolution', "data", "revis", "marker_stats_fine.Rdata"))
                    
# slurmjobs::job_single('revis01_find_markers_fine', create_shell = TRUE, memory = '25G', command = "Rscript revis01_find_markers_fine.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
