
library("RColorBrewer")

mdd_Dx_colors <- brewer.pal(3, "Dark2")
mdd_Dx_colors <- c("Control" = mdd_Dx_colors[[1]], #Green
                   "Bipolar" = mdd_Dx_colors[[2]], #Orange
                   "MDD" = mdd_Dx_colors[[3]]) #Purple
mdd_Dx_colors
# Control   Bipolar       MDD 
# "#1B9E77" "#D95F02" "#7570B3" 

mdd_Dx_colors_LD <- c("Control_light" = "#2EDCA7",
                      "Control_dark" = "#188C69",
                      "Bipolar_light" = "#FD9749",
                      "Bipolar_dark" = "#B65002",
                      "MDD_light" = "#9794C7",
                      "MDD_dark" = "#5A54A0")


mdd_BrainRegion_colors <- brewer.pal(3,"YlGnBu")

mdd_BrainRegion_colors <- c("Amygdala" = mdd_BrainRegion_colors[[1]],
                            "Overall" = mdd_BrainRegion_colors[[2]],
                            "sACC" = mdd_BrainRegion_colors[[3]])

mdd_dataset_colors <- c("psychENCODE_BP" = mdd_Dx_colors[[2]],
                        "psychENCODE_MDD" = mdd_Dx_colors[[3]])

# mdd_Sex_colors <- c(M = "#1E4482", `F` = "#783F58")
mdd_Sex_colors <- c(M = "steelblue4", `F` = "maroon4")

save(mdd_Dx_colors, mdd_BrainRegion_colors, mdd_Sex_colors, mdd_dataset_colors, file = "MDD_colors.Rdata")
