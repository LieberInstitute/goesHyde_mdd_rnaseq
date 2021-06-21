library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(here)

load(here("deconvolution","data","sce_filtered.Rdata"), verbose = TRUE)
ct <- list(split = "cellType.split", broad = "cellType")

cells_broad <- c("Inhib", "Excit","Astro","Micro","Oligo","OPC")
cells_inhib_sub <- c("Inhib.1","Inhib.2","Inhib.3","Inhib.4","Inhib.5")
cells_excit_sub <- c("Excit.1","Excit.2","Excit.3", "Excit.4")
cells_multi <- c("Multi-Excit", "Multi-Inhib", "Multi-Excit-Inhib")

cells_broad_colors <- c(brewer.pal(n = length(cells_broad), name = "Set1"))
names(cells_broad_colors) <- c(cells_broad)

inhib_colors <- colorRampPalette(c(cell_colors[["Inhib"]],"white"))(length(cells_inhib_sub)+1)
inhib_colors <- inhib_colors[1:length(cells_inhib_sub)]
names(inhib_colors) <- cells_inhib_sub

excit_colors <- colorRampPalette(c(cell_colors[["Excit"]],"white"))(length(cells_excit_sub)+1)
excit_colors <- excit_colors[1:length(cells_excit_sub)]
names(excit_colors) <- cells_excit_sub

multi_colors <- c("#1a3d59","#752262","#8a0f10")
names(multi_colors) <- cells_multi

multi_data <- tibble(cell_type_broad = c("Excit","Excit-Inhib","Inhib"), 
                     type = "Multi", cell_type = cells_multi)

cell_colors <- c(cells_broad_colors, inhib_colors, excit_colors,multi_colors)

test_data <- tibble(cell_type = c(cells_broad, cells_inhib_sub, cells_excit_sub)) %>%
  mutate(cell_type_broad = cell_type) %>%
  separate(cell_type_broad, "\\.", into = c("cell_type_broad", "type")) %>%
  mutate(type= ifelse(is.na(type), "Broad", as.character(type))) %>%
  rbind(multi_data)

test_data$type <- factor(test_data$type, levels = c("Broad","Multi","1","2","3","4","5"))

ggplot(test_data, aes(x= type, y= cell_type_broad, fill= cell_type))+
  geom_tile() +
  scale_fill_manual(values = cell_colors) +
  theme(legend.position = "None")

cell_colors
# Inhib             Excit             Astro             Micro             Oligo               OPC           Inhib.1           Inhib.2 
# "#E41A1C"         "#377EB8"         "#4DAF4A"         "#984EA3"         "#FF7F00"         "#FFFF33"         "#E41A1C"         "#E94749" 
# Inhib.3           Inhib.4           Inhib.5           Excit.1           Excit.2           Excit.3           Excit.4       Multi-Excit 
# "#EE7576"         "#F4A3A4"         "#F9D1D1"         "#377EB8"         "#699EC9"         "#9BBEDB"         "#CDDEED"         "#1a3d59" 
# Multi-Inhib Multi-Excit-Inhib 
# "#752262"         "#8a0f10"

save(cell_colors, file = "cell_colors.Rdata")
