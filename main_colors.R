
library("RColorBrewer")

mdd_Dx_colors <- brewer.pal(3, "Dark2")
mdd_Dx_colors <- c("Control" = mdd_Dx_colors[[1]],
                   "Bipolar" = mdd_Dx_colors[[2]],
                   "MDD" = mdd_Dx_colors[[3]])

mdd_BrainRegion_colors <- brewer.pal(3,"YlGnBu")

mdd_BrainRegion_colors <- c("Amygdala" = mdd_BrainRegion_colors[[1]],
                            "Overall" = mdd_BrainRegion_colors[[2]],
                            "sACC" = mdd_BrainRegion_colors[[3]])

mdd_dataset_colors <- c("psychENCODE_BP" = mdd_Dx_colors[[2]],
                        "psychENCODE_MDD" = mdd_Dx_colors[[3]])

mdd_Sex_colors <- c(M = "#1E4482", `F` = "#783F58")

mdd_all_colors <- c(mdd_Dx_colors, mdd_BrainRegion_colors, mdd_Sex_colors, mdd_dataset_colors)
