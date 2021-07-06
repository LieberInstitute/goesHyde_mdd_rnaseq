library(DeconvoBuddies)
library(here)

load(here("deconvolution","data","sce_filtered.Rdata"), verbose = TRUE)
(ct <- levels(sce_pan$cellType.Broad))

cell_colors <- create_cell_colors(pallet = "classic", cell_types = ct, preview = TRUE)
cell_colors
# Astro      Endo     Micro     Mural     Oligo       OPC     Tcell     Excit     Inhib
# "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628" "#F781BF" "#999999" "#E41A1C" "#377EB8"

cell_colors2 <- cell_colors
names(cell_colors2) <- c("Astro", "Micro", "Oligo", "OPC", "Mural","Endo", "Tcell", "Inhib", "Excit")
cell_colors <- cell_colors2[ct]
# Astro      Endo     Micro     Mural     Oligo       OPC     Tcell     Excit     Inhib
# "#4DAF4A" "#F781BF" "#984EA3" "#A65628" "#FF7F00" "#FFFF33" "#999999" "#377EB8" "#E41A1C"

save(cell_colors, file = here("deconvolution", "data","cell_colors.Rdata"))
