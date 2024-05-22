
library("DeconvoBuddies")
library("SingleCellExperiment")
library("purrr")
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
marker_stats <- map(sce_region, ~get_mean_ratio(.x, 
                               cellType_col = "cellType",
                               gene_ensembl = "gene_id",
                               gene_name = "gene_name")
)
                  
                    
                    
