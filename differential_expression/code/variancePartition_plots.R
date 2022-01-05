
library(tidyverse)
library(variancePartition)
library(sessioninfo)
library(here)

## Load and sort data
files <- list.files(path = here("differential_expression","data","variance_partition"), full.names = TRUE)
varPart <- lapply(files, function(x) mget(load(x, verbose = TRUE)))
varPart <- map(varPart, ~pluck(.x, 1))
names(varPart) <- gsub(".Rdata","",gsub("varPart_","",basename(files)))

## Sort columns
varPartSort <- map_depth(varPart, 2, sortCols)
map_depth(varPartSort, 2, colnames)

## Violin Plots
varPart_violin <- map2(varPartSort, names(varPartSort),  
                       function(data, data_name){
                         map2(data, names(data), ~plotVarPart(.x) +
                                labs(title = data_name,
                                     subtitle = paste("Model -",.y)) + 
                                theme(axis.text.x = element_text(angle = 90)))
                       })

pdf(here("differential_expression","plots", "variancePartition.pdf"))
map_depth(varPart_violin,2, print)
dev.off()
# 
# ## Density Plots
# varPartLong <- map(varPartSort, ~.x %>% 
#   as.data.frame() %>%
#   rownames_to_column("Gene") %>% 
#   pivot_longer(!Gene, names_to = "Term", values_to = "Variance"))
# 
# varPart_density <- map(varPartLong, ~ggplot(data = .x, aes(x = Variance, color = Term)) +
#                          geom_density() +
#                          )
# 
# map2(varPart_density, names(varPart_density), 
#      ~ggsave(.x, 
#              filename = here("differential_expression","plots",
#                              paste0("variancePartition_density_",.y,".png")), 
#              width = 10)
# )

## Compare Mean Variance 
varPart_t <- transpose(varPart)
names(varPart_t)

meanVar <- map(varPart_t, function(vp){
  
  meanVar <- do.call("cbind", map(vp, ~colMeans(.x))) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    pivot_longer(!term, values_to = "meanVar", names_to = "DE_analysis")
  
  term_order <- meanVar %>% 
    group_by(term) %>% 
    summarize(max = max(meanVar)) %>%
    arrange(-max) %>%
    pull(term)
  
  meanVar$term <- factor(meanVar$term, levels = term_order)
  return(meanVar)
})

## tile plots 
walk2(meanVar, names(meanVar), function(mV, name){
  
  meanVar_tile <- mV %>% 
    filter(term != "Residuals") %>%
    ggplot(aes(x = DE_analysis , y = term, fill = meanVar)) +
    geom_tile()  +
    scale_fill_continuous(type = "viridis") +
    labs(title = paste("Model -", name))
  
  ggsave(meanVar_tile, filename = here("differential_expression","plots", paste0("meanVar_tile_",name,".png")))
  
})

map(meanVar, ~.x %>% filter(term == "Residuals"))

