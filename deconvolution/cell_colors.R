library(RColorBrewer)

## load and melt prop data
load("prop_sacc.Rdata", verbose = TRUE)
load("prop_amyg.Rdata", verbose = TRUE)
load("prop_broad_sacc.Rdata", verbose = TRUE)
load("prop_broad_amyg.Rdata", verbose = TRUE)

## extract cell_types establish color pallet
cells_sacc <- colnames(est_prop_sacc$Est.prop.weighted)
cells_amyg <- colnames(est_prop_amyg$Est.prop.weighted)
cells_broad_sacc <- colnames(est_prop_broad_sacc$Est.prop.weighted)
cells_broad_amyg <- colnames(est_prop_broad_amyg$Est.prop.weighted)

cells_broad <- unique(c(cells_broad_amyg, cells_broad_sacc))
cells_specific <- unique(c(cells_sacc, cells_amyg))
cells_specific_unique <- cells_specific[!cells_specific %in% cells_broad]

cell_colors <- c(brewer.pal(n = length(cells_broad), name = "Set1"),
                 brewer.pal(n= length(cells_specific_unique), name = "Set3"))
names(cell_colors) <- c(cells_broad, cells_specific_unique)

cell_colors
# Inhib     Oligo     Astro     Excit     Micro       OPC   Inhib.2   Inhib.1   Excit.3   Excit.1   Excit.2   Excit.4 
# "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#8DD3C7" "#FFFFB3" "#BEBADA" "#FB8072" "#80B1D3" "#FDB462" 
# Inhib.3   Inhib.5   Inhib.4 
# "#B3DE69" "#FCCDE5" "#D9D9D9" 

rm(est_prop_sacc,est_prop_amyg, est_prop_broad_sacc, est_prop_broad_amyg)
