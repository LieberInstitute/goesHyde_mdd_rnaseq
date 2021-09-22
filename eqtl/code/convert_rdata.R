library(SummarizedExperiment)
library(sessioninfo)
library(tidyverse)
library(VariantAnnotation)
library(jaffelab)
library(here)

source(here("eqtl","code","rse_to_bed.R"))

## load
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_exon.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_jxn.Rdata"), verbose = TRUE)
load(here("exprs_cutoff", "rse_tx.Rdata"), verbose = TRUE)

## Split gene data
regions <- c(amyg = "Amygdala", sacc = "sACC")
features <- c("gene", "exon", "jxn", "tx")
names(features)  <- features

rse_gene_split <- map(regions, ~rse_gene[,rse_gene$BrainRegion == .x])
samples_split <- map(rse_gene_split, colnames)

#### Covariate Data ####
load(here("eqtl", "data", "pcs_4features.Rdata"), verbose = TRUE)
load(here("eqtl", "data", "pcs_4features_gene_temp.rda"), verbose = TRUE)
corner(genePCs)
dim(genePCs)

all(rownames(genePCs) == colnames(rse_gene))

covar_format <- function(data, rn){
  data <- as.data.frame(data) 
  rownames(data) <- rn
  data <- t(data)
  data <- as.data.frame(data) %>% rownames_to_column("id")
  return(data)
}

# covars <- map2(list(genePCs, exonPCs, jxnPCs, txPCs), features,
covars <- map2(list(gene = genePCs),features[1],
                function(pc, feat){
                  map2(rse_gene_split, regions, function(rse, region){
                    message(paste(feat, region))
                    ## Phenodata
                    pd <- colData(rse)[,c("PrimaryDx", "Sex", paste0("snpPC", 1:5))]
                    pd <- covar_format(pd, rse$genoSample)
                    
                    ## PC data
                    pc <- pc[colnames(rse),]
                    pc <- covar_format(pc, rse$genoSample)
                    
                    ##bind and save
                    covars <- rbind(pd, pc)
                    write.table(covars,
                                file= here("eqtl","data","covariates_txt", paste0("covariates_",feat,"_",region,".txt")),
                                sep = "\t", quote = FALSE, row.names = FALSE)
                    return(covars)
                  })
                })
# corner(covars$gene$amyg)


#### Expression Data ####
expression_fn <- map(features, function(feat) map(regions, ~here("eqtl","data", "expression_bed", paste0(feat, "_",.x,".bed"))))

expression_bed <- map2(list(rse_gene, rse_exon, rse_jxn, rse_tx), features, function(rse, feat){
  rse_split <- map(regions, ~rse[,rse$BrainRegion == .x])
  expr_bed <- map(rse_split, rse_to_bed)
  return(expr_bed)
})

walk2(expression_bed, expression_fn, function(expr, fn){
  walk2(expr, fn, ~write.table(.x, .y, sep = "\t", 
                               quote = FALSE, row.names = FALSE))
}
)

## Create shell script to zip data
commands <- map(unlist(expression_fn), ~paste0("bgzip ", .x," && tabix -p bed ", .x,".gz"))
if(file.exists("bed_bgzip.sh")) file.remove("bed_bgzip.sh")
sgejobs::job_single("bed_bgzip", create_shell = TRUE, queue= 'bluejay', memory = '100G', 
                    command = paste(commands, collapse = '\n'))
## Add "module load htslib"

#### VCF ####
risk_vcf <- readVcf(here("eqtl","data","risk_snps","LIBD_maf01_gwas_BPD.vcf.gz"))
risk_vcf
risk_vcf_split <- map(rse_gene_split, ~risk_vcf[,.x$genoSample])
map(risk_vcf_split, dim)

vcf_fn <- map(regions, ~here("eqtl","data","risk_snps",paste0("LIBD_maf01_gwas_BPD_",.x,".vcf.gz")))
walk2(risk_vcf_split, vcf_fn, ~writeVcf(.x, .y))

## plink commands
map(vcf_fn, ~paste("plink --make-bed --output-chr chrM --vcf", .x, "--out", gsub(".vcf.gz","",.x)))


## check 

map2(bed, risk_vcf_split, ~all(colnames(.x[,5:ncol(.x)]) == colnames(.y)))
map2(bed, covars, ~all(colnames(.x[,5:ncol(.x)]) == colnames(.y[[1]][,2:ncol(.y[[1]])])))

vcf_test <- readVcf("../data/risk_snps/LIBD_maf01_gwas_BPD_Amygdala.vcf.gz")

## prep interaction csv
walk2(rse_gene_split, regions, function(rse, region){
  cell_fractions <- colData(rse)[,c("Astro","Endo","Macro","Micro","Mural","Oligo","OPC","Tcell","Excit","Inhib")]
  cell_fractions <- as.data.frame(cell_fractions)
  rownames(cell_fractions) <- rse$genoSample
  write.csv(cell_fractions, file = here("eqtl","data","interaction",paste0("cell_fraction_",region,".csv")))
})

## create shell commands ##
sgejobs::job_single("tensorqtl_risk_snps", create_shell = TRUE, queue= 'bluejay', memory = '50G', 
                    command = "python tensorqtl_risk_snps.py")
