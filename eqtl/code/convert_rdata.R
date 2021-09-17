library(SummarizedExperiment)
library(sessioninfo)
library(tidyverse)
library(VariantAnnotation)
library(jaffelab)
library(here)

source(here("eqtl","code","rse_to_bed.R"))

## load
load(here("exprs_cutoff", "rse_gene.Rdata"), verbose = TRUE)
# load(here("exprs_cutoff", "rse_exon.Rdata"), verbose = TRUE)
# load(here("exprs_cutoff", "rse_jxn.Rdata"), verbose = TRUE)
# load(here("exprs_cutoff", "rse_tx.Rdata"), verbose = TRUE)

## Split gene data
regions <- c(amyg = "Amygdala", sacc = "sACC")
rse_gene_split <- map(regions, ~rse_gene[,rse_gene$BrainRegion == .x])

#### Covariate Data ####
# load(here("eqtl", "genomewide", "rdas", "pcs_4features_amyg.rda"), verbose = TRUE)
load(here("eqtl", "data", "pcs_4features_gene_temp.rda"), verbose = TRUE)

corner(genePCs)

all(rownames(genePCs) == colnames(rse_gene))

# map(list(genePCs, exonPCs, jxnPCs, txPCs), function(pc){
covars <- map2(rse_gene_split, regions,  ~map(list(genePCs), function(pc){
  message(.y)
  pc <- pc[colnames(.x),]
  rownames(pc) <- .x$genoSample
  pc <- t(pc)
  pc <- as.data.frame(pc) %>% rownames_to_column("id")
  write.table(pc,
              file= here("eqtl","data","covariates_txt", paste0("covariates_gene_",.y,".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  return(pc)
}))

ss(rownames(genePCs),"_")[(ss(rownames(genePCs),"_") %in% rse_gene$RNum)]


#### Create Phenotype Data ####
fn_gene <- map(regions, ~here("eqtl","data", "expression_bed", paste0("gene_",.x,".bed")))

bed <- map2(rse_gene_split, fn_gene, function(rse,fn){
  bed_gene <- rse_to_bed(rse)
  write.table(bed_gene, fn, sep = "\t", 
              quote = FALSE, row.names = FALSE)
  return(bed_gene)
})

## Create shell script to zip data
commands <- map(fn_gene, ~paste0("bgzip ", .x," && tabix -p bed ", .x,".gz"))
if(file.exists("bed_bgzip.sh")) file.remove("bed_bgzip.sh")
sgejobs::job_single("bed_bgzip", create_shell = TRUE, queue= 'bluejay', memory = '100G', 
                    command = paste(commands, collapse = '\n'))
## Add "module load htslib"

#### VCF ####
risk_vcf <- readVcf(here("eqtl","data","risk_snps","LIBD_maf01_gwas_BPD.vcf.gz"))
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
## create shell commands ##
