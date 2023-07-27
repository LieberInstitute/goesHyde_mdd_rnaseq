library("SummarizedExperiment")
library("sessioninfo")
library("tidyverse")
library("VariantAnnotation")
library("jaffelab")
library("here")
library("recount")

source(here("eqtl", "code", "rse_to_bed.R"))

## load
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_exon.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_jxn.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_tx.Rdata"), verbose = TRUE)

## add logcounts
assays(rse_gene)$logcounts <- log2(getRPKM(rse_gene, "Length")+1)
assays(rse_exon)$logcounts <- log2(getRPKM(rse_exon, "Length")+1)
assays(rse_jxn)$logcounts <- log2(getRPKM(rse_jxn, "Length")+1)
assays(rse_tx)$logcounts <- log2(assays(rse_tx)$tpm+1)

## Split gene data
regions <- c(amyg = "Amygdala", sacc = "sACC")
features <- c("gene", "exon", "jxn", "tx")
names(features) <- features

rse_gene_split <- map(regions, ~ rse_gene[, rse_gene$BrainRegion == .x])
samples_split <- map(rse_gene_split, colnames)

#### Covariate Data ####
pcs <- map(features, function(f_name) map(regions, ~ get(load(here("eqtl", "data", "featuresPCs", paste0(f_name, "PCs", .x, ".rda"))))))
corner(pcs$gene$amyg)

covar_format <- function(data, rn) {
    data <- as.data.frame(data)
    rownames(data) <- rn
    data <- t(data)
    data <- as.data.frame(data) %>% rownames_to_column("id")
    return(data)
}

covars <- map2(
    pcs, features,
    function(pc, feat) {
        pmap(list(rse = rse_gene_split, region = regions, pc = pc), function(rse, region, pc) {
            message(paste(feat, region))
            ## Phenodata
            pd <- as.data.frame(colData(rse)[, c("PrimaryDx", "Sex", paste0("snpPC", 1:5))])
            pd <- model.matrix(~ PrimaryDx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, data = pd)[, 2:9]
            pd <- covar_format(pd, rse$genoSample)

            ## PC data
            pc <- covar_format(pc, rse$genoSample)
            ## bind and save
            covars <- rbind(pd, pc)
            write.table(covars,
                file = here("eqtl", "data", "tensorQT_input", "covariates_txt", paste0("covariates_", feat, "_", region, ".txt")),
                sep = "\t", quote = FALSE, row.names = FALSE
            )
            return(covars)
        })
    }
)
corner(covars$gene$amyg)


#### Expression Data ####
expression_fn <- map(features, function(feat) map(regions, ~ here("eqtl", "data", "tensorQT_input", "expression_bed", paste0(feat, "_", .x, ".bed"))))

expression_bed <- map2(list(rse_gene, rse_exon, rse_jxn, rse_tx), features, function(rse, feat) {
    rse_split <- map(regions, ~ rse[, rse$BrainRegion == .x])
    expr_bed <- map(rse_split, rse_to_bed)
    return(expr_bed)
})

walk2(expression_bed, expression_fn, function(expr, fn) {
    walk2(expr, fn, ~ write.table(.x, .y,
        sep = "\t",
        quote = FALSE, row.names = FALSE
    ))
})

## Create shell script to zip data
commands <- map(unlist(expression_fn), ~ paste0("bgzip ", .x, " && tabix -p bed ", .x, ".gz"))
if (file.exists("02_bed_bgzip.sh")) file.remove("02_bed_bgzip.sh")
sgejobs::job_single("bed_bgzip",
    create_shell = TRUE, memory = "100G",
    command = paste(commands, collapse = "\n")
)
## Add "module load htslib"

#### VCF ####
risk_vcf <- readVcf(here("eqtl", "data", "risk_snps", "LIBD_maf01_gwas_BPD.vcf.gz"))
risk_vcf
risk_vcf_split <- map(rse_gene_split, ~ risk_vcf[, .x$genoSample])
map(risk_vcf_split, dim)

vcf_fn <- map(regions, ~ here("eqtl", "data", "risk_snps", paste0("LIBD_maf01_gwas_BPD_", .x, ".vcf.gz")))
walk2(risk_vcf_split, vcf_fn, ~ writeVcf(.x, .y))

## plink commands
map(vcf_fn, ~ paste("plink --make-bed --output-chr chrM --vcf", .x, "--out", gsub(".vcf.gz", "", .x)))


## check

# map2(bed, risk_vcf_split, ~ all(colnames(.x[, 5:ncol(.x)]) == colnames(.y)))
# map2(bed, covars, ~ all(colnames(.x[, 5:ncol(.x)]) == colnames(.y[[1]][, 2:ncol(.y[[1]])])))


## prep interaction csv
walk2(rse_gene_split, regions, function(rse, region) {
    cell_fractions <- colData(rse)[, c("Astro", "Endo", "Macro", "Micro", "Mural", "Oligo", "OPC", "Tcell", "Excit", "Inhib")]
    cell_fractions <- as.data.frame(cell_fractions)
    rownames(cell_fractions) <- rse$genoSample
    write.csv(cell_fractions, file = here("eqtl", "data", "tensorQT_input", "interaction", paste0("cell_fraction_", region, ".csv")))
})

## create shell commands ##
# sgejobs::job_single("tensorqtl_risk_snps",
#     create_shell = TRUE, queue = "bluejay", memory = "50G",
#     command = "python tensorqtl_risk_snps.py"
# )

# sgejobs::job_single('01_convert_rdata_to_tensorqtl_format', create_shell = TRUE, memory = '50G', command = "Rscript 01_convert_rdata_to_tensorqtl_format.R")

## Reproducibility information
print("Reproducibility information:")
options(width = 120)
session_info()
