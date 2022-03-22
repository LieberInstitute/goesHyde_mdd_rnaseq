eQTL Analysis
========
Using tensorQTL https://github.com/broadinstitute/tensorqtl

## Risk SNP interactions 
Using `cis.map_nominal()`  
Only ran Gene level  
Important p-value: "**pval_gi** - p-value of gi in the model" (tensorqtl README)  

### BPD Risk SNPs
VCF: `./data/risk_snps/LIBD_maf01_gwas_BDP.vcf.gz` was made by Josh Stolz, 9/15/2021  
Contains 63/64 risk SNPs from: `./data/risk_snps/pgc3-supptable2.csv` (Table2 from pgc3 BipSeq supplementary tables)  
tensorQTL code: `./code/tensorqtl_risk_interaction_BPD.py`  
Raw Output: `./data/tensorQTL_out/nominal_bpd_risk/`  csv file for each cell type
Summary Output: `./data/summary/gene_BPD_risk_cell_fraction_interaction.csv` & `./data/summary/gene_BPD_risk_cell_fraction_interaction_summary.csv`  

### MDD Risk SNPs
VCF: `./data/risk_snps/LIBD_maf01_gwas_MDD.vcf.gz`  
Contains 2152/4625 SNPs from:`./data/risk_snps/PGC_depression_genome-wide_significant_makers.txt`  
tensorQTL code: `./code/tensorqtl_risk_interaction_MDD.py`
Raw Output: `./data/tensorQTL_out/nominal_mdd_risk/`  csv file for each cell type
Summary Output: `./data/summary/gene_MDD_risk_cell_fraction_interaction.csv` & `./data/summary/gene_MDD_risk_cell_fraction_interaction_summary.csv`  


## Genomewide
VCF:  `../genotype_data/mdd_bpd/maf01/mdd_bpd_maf01.vcf.gz` ~11M SNPs

### Nominal - all results
Using `cis.map_nominal()`  
Ran All 4 Features  
Important p-val: "**pval_nominal** - nominal p-value of association" (FASTQTL docs)  
tensorQTL code: `./code/tensorqtl_genomewide_nominal_argv.py`  
Raw Output: `./data/tensorQTL_out/genomewide_nominal`  
  * Contains parquet files   
Filtering code to make FDR csv files: `./code/filter_genomewide_nominal.R`
  * `.data/tensorQTL_FDR01/geomewide_nominal/`

Risk SNP subsets: `./data/risk_snps_eqtl/` One csv for each risk SNP set 
(MDD or BPD) + region + feature.  
Code to subset: `./code/subset_genomewide_nominal_riskSNPs.R`
Summary code: `./code/summarize_genomewide_nominal.R`  

### Cis - best SNP per feature (gene)  
Using `cis.map_cis()`  
Ran Gene Only  
Important p-val: "**pval_beta** - permutation p-value obtained via beta approximation." (FASTQTL docs)   
We advice to use this one in any downstream analysis."  
qval column: FDR corrected `pval_beta` added for map_independent run
tensorQTl code: `./code/tensorqtl_genomewide_cis.py`  
Summary code: `./code/summarize_genomewide_cis.R`  

### Independent - conditionally independent cis-QTLs using the stepwise regression procedure  
Using Run `cis.map_independent()`  
Ran Gene Only  
Window = 500 KB  
Cis results are input, FDR cutoff 0.01  
tensorQTl code: `./code/tensorqtl_genomewide_independent.sh`
runs `./code/tensorqtl_genomewide_independent.py`  on the GPU!  
Summary code: `./code/summarize_genomewide_independent.R`  


## DataPrep  
`./code/convert_rdata.R`: convert data to tables needed for tensorQTL  
* also generates plink commands if needed

## Input Data
Covariate data (from the `colData`): `./data/tensorQTL_input/covariates_txt`  
Expression data (from the `colData`): `./data/tensorQTL_input/expression_bed`  
Interaction data (cell type proportion estimates): `./data/tensorQTL_input/covariates_txt`  

## Output Data
Raw output from tenosrQTL: `./data/tensorQTL_out`  
Filtered output for FDR < 0.01: `./data/tensorQTL_FDR01`  

## NOTES:  
Some reorganization may cause filenames to not match up perfectly...  
`genomewide/` is an old folder containing matrixQTL scripts, maybe remove someday
