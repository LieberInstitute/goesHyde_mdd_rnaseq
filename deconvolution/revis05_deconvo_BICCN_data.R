
library("SingleCellExperiment")
library("spatialLIBD")
library("Seurat")
library("hspe")
library("tidyverse")
library("DeconvoBuddies")
library("here")
library("sessioninfo")

## get args
args = commandArgs(trailingOnly=TRUE)
region <- args[1]
# region = "Amygdala"

#### load data
message(Sys.time(), " - Load rse, filter to ", region)

load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
rownames(rse_gene) <- rowData(rse_gene)$ensemblID

table(rse_gene$BrainRegion)

rse_gene <- rse_gene[,rse_gene$BrainRegion == region]
ncol(rse_gene)

message(Sys.time(), " - Load sce object")
# Dissection: Amygdaloid complex (AMY) - Central nuclear group - CEN
# Dissection: Cerebral cortex (Cx) - Subcallosal Gyrus (SCG) - Subgenual (subcallosal) division of MFC - A25
sce_path <- c(Amygdala = here("deconvolution", "BICCN_data", "AMY_CEN", "c821776d-a9fd-4a9b-b0ef-65bee87c7dbc.rds"),
                 sACC = here("deconvolution", "BICCN_data", "SCG_MFC", "82faf671-658a-4c88-8a7d-618fc7f68fad.rds"))[region]

print(sce_path)

sce <- Seurat::as.SingleCellExperiment(readRDS(sce_path))
dim(sce)
# 59236 42699

table(sce$donor_id)

(n_ct <- table(sce$supercluster_term))

## Drop small (<100) and ambiguous cell types 
drop_ct <- c(names(n_ct)[n_ct < 100], "Splatter", "Miscellaneous")
message("dropping: ", paste0(drop_ct, collapse = ", ") )
sce <- sce[, !sce$supercluster_term %in% drop_ct]
sce$supercluster_term <- droplevels(sce$supercluster_term)

## make cell type names syntactically valid
sce$supercluster <- factor(gsub(" |-", "_", sce$supercluster_term))

## save cell type proportions 

sce_cell_type_prop <- as.data.frame(colData(sce)) |>
  count(donor_id, supercluster) |>
  group_by(donor_id) |>
  mutate(prop = n/sum(n))

write.csv(sce_cell_type_prop, file = here("deconvolution", "BICCN_data", paste0("cell_type_prop_", region, ".csv")), row.names = FALSE)

## filter to common genes
common_genes <- intersect(rownames(rse_gene), rownames(sce))
sce <- sce[common_genes, ]

message("\nSce dim:")
dim(sce)

table(sce$supercluster)

#### Find Marker Genes ####
message(Sys.time(), " - Find marker genes")
marker_stats <- get_mean_ratio(sce, cellType_col = "supercluster")

save(marker_stats, file = here("deconvolution", "BICCN_data", paste0("marker_stats_", region, ".Rdata")))

#### pseudbulk sce data for hspe ####
message(Sys.time(), " - pseudobulk data")

## logcounts looks like the counts... (all ints)
counts(sce) <- logcounts(sce)
logcounts(sce) <- NULL

sce_pb <- spatialLIBD::registration_pseudobulk(
  sce,
  var_registration = "supercluster",
  var_sample_id = "donor_id",
  covars = NULL,
  min_ncells = 10,
  pseudobulk_rds_file = NULL
)
dim(sce_pb)

rm(sce)

#### Run hspe ####
message(Sys.time(), " - Prep hspe input")
pure_samples = rafalib::splitit(sce_pb$supercluster)

common_genes <- intersect(rownames(rse_gene), rownames(sce_pb))
message("common genes: ", length(common_genes))

## extract count data
mixture_samples = t(log2(assays(rse_gene)$counts+1)[common_genes,])
reference_samples = t(assays(sce_pb)$logcounts[common_genes,])

stopifnot(ncol(mixture_samples) == ncol(reference_samples))

## filter marker stats
marker_stats <- marker_stats |>
  filter(MeanRatio.rank <= 25 & MeanRatio > 1 & gene %in% common_genes)

write.csv(marker_stats, file = here("deconvolution", "BICCN_data", paste0("marker_stats_top25_", region, ".csv")), row.names = FALSE)

marker_stats |> count(cellType.target)

marker_genes <- purrr::map(rafalib::splitit(marker_stats$cellType.target), ~marker_stats$gene[.x])
marker_genes <- marker_genes[names(pure_samples)]

message(Sys.time(), " - run hspe")
est_prop_hspe = hspe(Y = mixture_samples,
                reference = reference_samples,
                pure_samples = pure_samples,
                markers = marker_genes,
                seed = 10524)


## Add long data and save
message(Sys.time(), " - Done! Saving!")

## save
save(est_prop_hspe, file = here("deconvolution", "BICCN_data",paste0("est_prop_BICCN_hspe_",region,".Rdata")))

# slurmjobs::job_loop(loops = list(region = c("Amygdala", "sACC")),
#                     name = 'revis05_deconvo_BICCN_data', 
#                     create_shell = TRUE, 
#                     memory = '25G')

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

