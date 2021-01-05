library(SummarizedExperiment)
library(RColorBrewer)
library(jaffelab)
library(here)
library(reshape2)
library(patchwork)
library(purrr)
library(tidyverse)
library(broom)
library(sessioninfo)
library(pheatmap)
library(SingleCellExperiment)

## Load colors and plotting functions
source(here("main_colors.R"))
source(here("deconvolution","big_little_boxplot.R"))
load(here("deconvolution","data","cell_colors.Rdata"), verbose = TRUE)

## Load MuSiC results & res data
load(here("deconvolution","data","est_prop_top5_long.Rdata"), verbose = TRUE)
load(here("deconvolution","data","est_prop_top5.Rdata"),verbose = TRUE)
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)

#### Boxplots ####
walk2(est_prop_long, names(est_prop_long), function(x,y){
  blb <- big_little_boxplot(data = x,
                            xvar = "cell_type", 
                            yvar = "prop",
                            fillvar =  "PrimaryDx",
                            colorvar = "ignore",
                            title = "MuSiC Proptions: Top5 markers",
                            subtitle = y)
  ggsave(plot = blb, filename = paste0("plots/cellType_boxplots_",y,".png"), width = 15)
} )


dx_sex_colors <- mdd_Dx_colors_LD
names(dx_sex_colors) <- gsub("dark","F",gsub("light","M",names(mdd_Dx_colors_LD)))

walk2(est_prop_long, names(est_prop_long), function(x,y){
  
  x <- x %>% mutate(Dx_Sex = paste0(PrimaryDx, "_",Sex))
  
  blb <- big_little_boxplot(data = x,
                            xvar = "cell_type", 
                            yvar = "prop",
                            fillvar =  "Dx_Sex",
                            colorvar = "ignore",
                            pallet = dx_sex_colors,
                            title = "MuSiC Proptions: Top5 markers",
                            subtitle = y)
  ggsave(plot = blb, filename = paste0("plots/cellType_sex_boxplots_",y,".png"), width = 15)
} )


#### Experimental plots ####
bar <-est_prop_long[[1]]%>%
  left_join(est_prop_long[[1]] %>%
              filter(cell_type == 'Oligo') %>%
              select(sample, prop_1 = prop),
            by = 'sample') %>% 
  ggplot(aes(x = reorder(sample, desc(prop_1)), y = prop, fill = cell_type)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cell_colors)+
  coord_flip()

ggsave(bar, filename = "plots/bar_test.png", height = 10)

## try pheatmap
my_anno_colors <- list(cellType = cell_colors,
                       PrimaryDx = mdd_Dx_colors,
                       Experiment = mdd_dataset_colors)

sample_anno <- as.data.frame(colData(rse_gene)) %>%
  select(Experiment, Sex, PrimaryDx, AgeDeath)

walk2(est_prop, names(est_prop), function(est, n){
  epw <- est$Est.prop.weighted
  png(paste0("plots/MuSiC-heatmap_",n,".png"), height = 800, width = 580)
  pheatmap(epw,
           show_rownames = FALSE,
           annotation_row = sample_anno,
           annotation_colors = my_anno_colors,
           main = paste(n,"cell type proptions"))
  dev.off()
})

#### Compare MuSiC output with sn data from same donors ####
load(here("deconvolution","data","sce.sacc_filtered.Rdata"), verbose = TRUE)
load(here("deconvolution","data","sce.amyg_filtered.Rdata"), verbose = TRUE)

donors <- unique(sce.sacc$donor)

input_cells <- list(sacc_broad = colData(sce.sacc)[,c("donor","cellType.Broad")],
                   amyg_broad = colData(sce.amyg)[,c("donor","cellType.Broad")],
                   sacc_specific = colData(sce.sacc)[,c("donor","cellType")],
                   amyg_specific = colData(sce.amyg)[,c("donor","cellType")])

input_prop <- map(input_cells, function(x){
  colnames(x)[2] <- "cell_type"
  
  x %>% as_tibble %>%
    group_by(donor, cell_type) %>%
    count() %>%
    group_by(donor) %>%
    mutate(input_prop = n/sum(n))
})

output_prop <- map2(est_prop_long, input_prop, function(est,input){
  est %>% select(donor = BrNum, cell_type, prop) %>%
    inner_join(input)
})

output_prop_df <- do.call("rbind", output_prop) %>%
  rownames_to_column(var = "dataset") %>%
  mutate(dataset = ss(dataset,"\\."))

in_out_prop <- ggplot(output_prop_df, aes(input_prop, prop, color = cell_type, shape = donor))+
  geom_point() +
  scale_color_manual(values = cell_colors)+
  geom_abline()+
  facet_wrap(~dataset) +
  labs(y = "MuSiC estimated prop", x = "Input sn data prop")

ggsave(in_out_prop, filename = "plots/prop_input_vs_ouput.png")

# sgejobs::job_single('deconvo_plots', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript deconvo_plots.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessionsession_info()
