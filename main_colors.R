
library("RColorBrewer")

mdd_colors <- brewer.pal(3, "Dark2")
mdd_colors <- c("Control" = mdd_colors[[1]],
               "Bipolar" = mdd_colors[[2]],
               "MDD" = mdd_colors[[3]])
