
library(SingleCellExperiment)
library(reshape2)
library(tidyverse)

get_mean_ratio <- function(sce, cellType_col =  "cellType"){
  
  sce_celltypes <- as.data.frame(colData(sce)) %>%
    select(cellType = !!sym(cellType_col)) %>%
    rownames_to_column(var = "id") 
  
   gene_stat <- as.matrix(assays(sce)$logcounts) %>%
     melt() %>%
     rename(gene = Var1, id = Var2, logcounts = value) %>%
     left_join(sce_celltypes, by = "id")  %>%
     group_by(gene, cellType) %>%
     summarise(median_logcount = median(logcounts),
               mean_logcount = mean(logcounts)) %>%
     ungroup()
  
   target_stat <- gene_stat %>% 
     filter(median_logcount != 0) %>%
     rename(cellType.target = cellType,
            mean_logcount.target = mean_logcount,
            median_logcount.target = median_logcount)
   
   mean_prop <- target_stat %>% right_join(gene_stat, by = "gene") %>%
     filter(cellType.target != cellType) %>%
     mutate(ratio = (mean_logcount.target + 0.01)/(mean_logcount + 0.01)) %>%
     arrange(gene, cellType.target,ratio) %>%
     group_by(gene, cellType.target) %>%
     slice(1) %>%
     group_by(cellType.target) %>%
     arrange(-ratio) %>%
     mutate(ratio_rank = row_number())
   
   return(mean_prop)
}
