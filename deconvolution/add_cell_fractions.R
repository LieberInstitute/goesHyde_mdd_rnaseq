library(SummarizedExperiment)
library(here)

load(here("deconvolution","data","est_prop_Bisque.Rdata"),verbose = TRUE)
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)

est_prop_bisque <- est_prop_bisque$all

all(colnames(rse_gene) == rownames(est_prop_bisque$bulk.props))

colData(rse_gene) <- cbind(colData(rse_gene), est_prop_bisque$bulk.props)
colnames(colData(rse_gene))

save(rse_gene, file = here("exprs_cutoff", "rse_gene.Rdata"))

# sgejobs::job_single('add_cell_fractions', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript add_cell_fractions.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
