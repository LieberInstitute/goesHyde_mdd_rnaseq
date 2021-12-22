
library(tidyverse)
library(variancePartition)
library(sessioninfo)
library(here)

## Load and sort data
files <- list.files(path = here("differential_expression","data","variance_partition"), full.names = TRUE)
# results <- sapply(files, function(x) mget(load(x)), simplify = TRUE) 

## for now
load(files[[4]], verbose = TRUE)
varPartSort <- map(varPart, sortCols)
map(varPartSort, colnames)

## Violin Plots
varPart_violin <- map2(varPartSort, names(varPartSort), ~plotVarPart(.x) + labs(title = paste("Model -",.y)))

map2(varPart_violin, names(varPart_violin), 
     ~ggsave(.x, 
             filename = here("differential_expression","plots",
                             paste0("variancePartition_",.y,".png")), 
             width = 10)
     )

pdf(here("differential_expression","plots", "variancePartition.pdf"))
walk(varPart_violin, print)
dev.off()

## Density Plots
varPartLong <- map(varPartSort, ~.x %>% 
  as.data.frame() %>%
  rownames_to_column("Gene") %>% 
  pivot_longer(!Gene, names_to = "Term", values_to = "Variance"))

varPart_density <- map(varPartLong, ~ggplot(data = .x, aes(x = Variance, color = Term)) +
                         geom_density() +
                         )

map2(varPart_density, names(varPart_density), 
     ~ggsave(.x, 
             filename = here("differential_expression","plots",
                             paste0("variancePartition_density_",.y,".png")), 
             width = 10)
)


