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
  * Contains parquet and FDR01 Summary csv files   
  * filtering code to make csv files: `./code/filter_genomewide_nominal.R`
Summary counts: `./data/summary/cis_genomewide_FDR01.csv`  
Risk SNP subsets: `./data/risk_snps_eqtl/` One csv for each risk SNP set 
(MDD or BPD) + region + feature. Code to subset: `./code/subset_genomewide_nominal_riskSNPs.R`   

### Cis - best SNP per feature
Using `cis.map_cis()`  
Ran Gene Only  
Important p-val: "**pval_beta** - permutation p-value obtained via beta approximation." (FASTQTL docs)   
We advice to use this one in any downstream analysis."  
qval column: FDR corrected `pval_beta` added for map_independent run
tensorQTl code: `./code/tensorqtl_genomewide_cis.py`  
Raw Output: `./data/tensorQTL_out/genomewide_cis`  
Summary counts: TODO  
TODO: Run `cis.map_independent()`

## DataPrep  
`./code/convert_rdata.R`: convert data to tables needed for tensorQTL  
* Covaraite data  
* Expression Data  
* plink commands if needed  
* Interaction aka cell fraction data  


## NOTES:  
Some reorganization may cause filenames to not match up perfectly...  
`genomewide/` is an old folder containing matrixQTL scripts, maybe remove someday
