library(SingleCellExperiment)
# library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
# library(batchelor)
# library(DropletUtils)
library(jaffelab)
library(limma)
library(here)
library(reshape2)
library(stringr)
library(purrr)
library(tidyverse)
#### Load and filter data ####
## Load rse_gene data
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

## sACC data
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_sACC-n2_cleaned-combined_SCE_MNTFeb2020.rda", verbose = TRUE)
sce.sacc$uniqueID <- paste0(sce.sacc$donor, "_", sce.sacc$Barcode)
colnames(sce.sacc) <- sce.sacc$uniqueID

sce.sacc <- sce.sacc[,sce.sacc$cellType != "Ambig.lowNtrxts",]
sce.sacc$cellType <- droplevels(sce.sacc$cellType)
## Add cellType.broad
sce.sacc$cellType.Broad <- ss(as.character(sce.sacc$cellType), "\\.", 1)
sce.sacc$cellType.Broad <- factor(sce.sacc$cellType.Broad)
## Match rownames
rownames(sce.sacc) <- rowData(sce.sacc)$ID
table(rownames(sce.sacc) %in% rownames(rse_gene))
# FALSE  TRUE 
# 15021 18517 

## filter for common between sc and bulk data
common_genes <- rownames(rse_gene)[rownames(rse_gene) %in% rownames(sce.sacc)]
sce.sacc <- sce.sacc[common_genes, ]

# Remove 0 genes across all nuclei
sce.sacc <- sce.sacc[!rowSums(assay(sce.sacc, "counts"))==0, ] 
dim(sce.sacc)
# [1] 17785  7004

#### Find markers ####
## create filters for median != 0
cellSubtype.idx <- splitit(sce.sacc$cellType.Broad)
medianNon0.idx <- lapply(cellSubtype.idx, function(x){
  apply(as.matrix(assay(sce.sacc, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

sapply(medianNon0.idx, sum)
# Astro Excit Inhib Micro Oligo   OPC 
# 1501  6892  5658   992  1259  2279

## Traditional t-test with design as in PB'd/limma approach ===
mod <- with(colData(sce.sacc), model.matrix(~ donor))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

markers.sacc.t.1vAll <- list()
for(i in levels(sce.sacc$cellType.Broad)){
  # Make temporary contrast
  message(i)
  sce.sacc$contrast <- ifelse(sce.sacc$cellType.Broad==i, 1, 0)
  # Test cluster vs. all
  markers.sacc.t.1vAll[[i]] <- findMarkers(sce.sacc, groups=sce.sacc$contrast,
                                           assay.type="logcounts", design=mod, test="t",
                                           direction="up", pval.type="all", full.stats=T)
}

temp.1vAll <- list()
for(i in levels(sce.sacc$cellType.Broad)){
  # Make temporary contrast
  message(i)
  sce.sacc$contrast <- ifelse(sce.sacc$cellType.Broad==i, 1, 0)
  # Test cluster vs. all
  temp.1vAll[[i]] <- findMarkers(sce.sacc, groups=sce.sacc$contrast,
                                 assay.type="logcounts", design=mod, test="t",
                                 std.lfc=TRUE,
                                 direction="up", pval.type="all", full.stats=T)
}

# Replace that empty slot with the entry with the actual stats
# What is stats 1 vs. stats 0? 
markers.sacc.t.1vAll <- lapply(markers.sacc.t.1vAll, function(x){ x[[2]] })
# Same for that with std.lfc
temp.1vAll <- lapply(temp.1vAll, function(x){ x[[2]] })

# Save all marker stats
m <- markers.sacc.t.1vAll
# Now just pull from the 'stats.0' DataFrame column
markers.sacc.t.1vAll <- lapply(markers.sacc.t.1vAll, function(x){ x$stats.0 })
temp.1vAll <- lapply(temp.1vAll, function(x){ x$stats.0 })

# Re-name std.lfc column and add to the first result
for(i in levels(sce.sacc$cellType.Broad)){
  colnames(temp.1vAll[[i]])[1] <- "std.logFC"
  markers.sacc.t.1vAll[[i]] <- cbind(markers.sacc.t.1vAll[[i]], temp.1vAll[[i]]$std.logFC, log10(m[[i]]$p.value))
  # Oh the colname is kept weird
  colnames(markers.sacc.t.1vAll[[i]])[4:5] <- c("std.logFC","log10.p.value")
  # Then re-organize
  markers.sacc.t.1vAll[[i]] <- markers.sacc.t.1vAll[[i]][ ,c("logFC","std.logFC","log.p.value","log10.p.value","log.FDR")]
  markers.sacc.t.1vAll[[i]] <- cbind(markers.sacc.t.1vAll[[i]],
                                     medianNon0.idx[[i]][match(rownames(markers.sacc.t.1vAll[[i]]),
                                                               names(medianNon0.idx[[i]]))])
  colnames(markers.sacc.t.1vAll[[i]])[6] <- "non0median"
}

save(markers.sacc.t.1vAll, file = here("deconvolution","data","markers_broad.sacc.Rdata"))

#### Classify and add stats ####

# load(here("deconvolution","data","markers_broad.sacc.Rdata"), verbose = TRUE)
markerList.t.1vAll <- lapply(markers.sacc.t.1vAll, function(x){
  rownames(x)[x$log.FDR < -6 & x$non0median==TRUE]
}
)

lengths(markerList.t.1vAll)
# Astro Excit Inhib Micro Oligo   OPC 
# 832  5162  4186   708   841  1160    
map_int(markers.sacc.t.1vAll, ~sum(.x$log10.p.value == -Inf))
# Astro Excit Inhib Micro Oligo   OPC 
# 249   789   226   322   360   148 
top_n = 40 
genes.top.t <- lapply(markerList.t.1vAll, function(x){head(x, n=top_n)})

## Save new markers
broad_markers_sacc <-data.frame(gene = unlist(markerList.t.1vAll))
broad_markers_sacc$cellType.Broad <- gsub("\\d+","",rownames(broad_markers_sacc))
write.csv(broad_markers_sacc,file = "data/braod_markers_sacc.csv")

## Find ratios
sce_celltypes <- as.data.frame(colData(sce.sacc)) %>%
  select(uniqueID, cellType.Broad)

marker_medians <- as.matrix(assays(sce.sacc)$logcounts[unlist(markerList.t.1vAll),]) %>%
  melt() %>%
  rename(gene = Var1, uniqueID = Var2, logcounts = value) %>%
  left_join(sce_celltypes, by = "uniqueID")  %>%
  group_by(gene,cellType.Broad) %>%
  summarise(median_log_count = median(logcounts))%>%
  arrange(gene, -median_log_count) %>%
  mutate(rank = row_number()) %>%
  left_join(broad_markers_sacc %>% rename(cellType.marker = cellType.Broad), by = "gene") %>%
  ungroup()

marker_medians %>% ungroup() %>% filter(cellType.Broad == cellType.marker) %>% count(rank)

marker_ratio <- marker_medians %>% filter(cellType.Broad == cellType.marker) %>%
  select(gene, cellType.marker, marker_median = median_log_count) %>%
  left_join(marker_medians %>% filter(cellType.Broad != cellType.marker)) %>%
  mutate(ratio = (marker_median + 0.01)/(median_log_count + 0.01),
         anno = paste0(cellType.marker,"/",cellType.Broad," = ",round(ratio, 3))) %>%
  group_by(gene, cellType.marker) %>%
  slice(1) %>%
  group_by(cellType.marker) %>%
  arrange(-ratio) %>%
  mutate(ratio_rank = row_number()) %>%
  ungroup() %>%
  select(-rank)

## compare with other stats
markers.sacc.table <- do.call("rbind",markers.sacc.t.1vAll) %>% 
  as.data.frame() %>%
  rownames_to_column("gene")%>%
  as_tibble() %>%
  add_column(cellType.marker = rep(names(markers.sacc.t.1vAll), each = nrow(sce.sacc))) %>%
  mutate(gene = ss(gene,"\\.")) %>%
  group_by(cellType.marker) %>%
  mutate(marker_rank = row_number())

markers.sacc.table$Symbol <- rowData(sce.sacc)[markers.sacc.table$gene,]$Symbol

marker_stat <- markers.sacc.table %>%
  right_join(marker_ratio, by = c("gene", "cellType.marker")) %>%
  ungroup() %>%
  mutate(Feature = paste0(str_pad(ratio_rank, 4, "left"),": ", gene, "-", Symbol),
         anno = paste0(" ",anno,"\n std logFC = ", round(std.logFC,3)))

#### Create Ratio plots #### 
ratio_plot <- ggplot(marker_stat, aes(x=ratio, y=std.logFC, color = cellType.Broad))+
  geom_point(size = .5) +
  facet_wrap(~cellType.marker, scales = "free_x") +
  labs(x = "target + 0.01/highest non-target + 0.01")+ 
  scale_x_continuous(trans = 'log10') +
  scale_color_manual(values = cell_colors)

ggsave(filename = "plots/expr/ratio_vs_stdFC.png")

ratio_median_plot <- ggplot(marker_stat, aes(x=ratio, y=marker_median))+
  geom_point(size = .5) +
  facet_wrap(~cellType.marker, scales = "free_x") +
  labs(x = "target + 0.01/highest non-target + 0.01")

ggsave(plot= ratio_median_plot, filename = "plots/expr/ratio_vs_median.png")

#### Plot expression ####
load("data/cell_colors.Rdata", verbose = TRUE)
for(i in names(markerList.t.1vAll)){
  marker_stat_cell <- marker_stat %>% filter(cellType.marker == i, ratio_rank <= 40)
  temp_sce <- sce.sacc[marker_stat_cell$gene,]
  rownames(temp_sce) <- marker_stat_cell$Feature
  
  png(paste0("plots/expr/sACC_t-sn-level_1vALL_topMarkers-REFINED-",i,"_logExprs.png"), height=1900, width=1200)
  print(
    plotExpression(temp_sce, exprs_values = "logcounts", features = marker_stat_cell$Feature,
                   x="cellType.Broad", colour_by="cellType.Broad", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + 
      scale_color_manual(values = cell_colors)+
      stat_summary(fun = median, fun.min = median, fun.max = median,
                   geom = "crossbar", width = 0.3) +
      geom_text(data = marker_stat_cell, aes(x = -Inf, y = Inf, label = anno),
                vjust = "inward", hjust = "inward")+
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      # ggtitle(label=paste(i, "top", top_n ,"markers, refined: single-nucleus-level p.w. t-tests, cluster-vs-all")
      ggtitle(label=paste(i, "top", top_n ,"markers: Ranked by median ratios"))
  )
  dev.off()
}

#### Check Overlaps with old data ####
top40_old_sacc <- read.csv("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/tables/top40genesLists_sACC-n2_cellType_SN-LEVEL-tests_May2020.csv")
top40_old_sacc <- top40_old_sacc[,grepl("1vAll", colnames(top40_old_sacc))]

top40_new_sacc <- read.csv("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/tables/top40genesLists-REFINED_sACC-n2_cellType.split_Nov2020.csv")
top40_new_sacc <- top40_new_sacc[,grepl("1vAll", colnames(top40_new_sacc))]

data.frame(cell_types = ss(colnames(top40_old_sacc),"_"),
           overlap_old = map_int(names(top40_old_sacc), ~sum(top40_old_sacc[[.x]] %in% rowData(sce.sacc)[genes.top.t[[ss(ss(.x,"_"),"\\.")]],]$Symbol)),
           overlap_new = map_int(names(top40_new_sacc), ~sum(top40_new_sacc[[.x]] %in% rowData(sce.sacc)[genes.top.t[[ss(ss(.x,"_"),"\\.")]],]$Symbol)))
#    cell_types overlap_old overlap_new
# 1       Astro          40          40
# 2     Excit.1          13          13
# 3     Excit.2           8           8
# 4     Excit.3           3           4
# 5     Excit.4           0           0
# 6     Inhib.1          20          20
# 7     Inhib.2          15          18
# 8       Micro          38          38
# 9       Oligo          39          39
# 10        OPC          37          37 


