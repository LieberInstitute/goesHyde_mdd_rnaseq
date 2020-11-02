library(ggplot2)
library(patchwork)
library("rlang")

big_little_boxplot <- function(data, xvar, yvar, fillvar, pallet = mdd_Dx_colors,  cutoff = 0.1, title, subtitle, fn){
  
  data <- data %>% group_by(cell_type) %>% mutate(big = mean(prop) > 0.1)
  
  bp_big <- data %>% filter(big) %>%
    ggplot(aes(x = !!sym(xvar), y = !!sym(yvar), fill = !!sym(fillvar))) +
    geom_boxplot()+
    scale_fill_manual(values = pallet)+
    labs(title = title,
         subtitle = subtitle)+
    theme(legend.position = "None")
  
  bp_little <- data %>% filter(!big) %>%
    ggplot(aes(x = !!sym(xvar), y = !!sym(yvar), fill = !!sym(fillvar))) +
    geom_boxplot()+
    scale_fill_manual(values = pallet)
  
  ggsave(filename = fn,plot = bp_big + bp_little, width = 10)
}
