# eQTL Analysis
Using tensorQTL https://github.com/broadinstitute/tensorqtl

## Risk SNP interactions 
Using `cis.map_nominal()`

### BPD Risk SNPs
VCF: `./data/risk_snps/LIBD_maf01_gwas_BDP.vcf.gz` was made by Josh Stolz, 9/15/2021  
Contains 63/64 risk SNPs from: `./data/risk_snps/pgc3-supptable2.csv` (Table2 from pgc3 BipSeq supplementary tables)  
tensorQTL code: `./code/tensorqtl_risk_interaction_BPD.py`  
Raw Output: `./data/tensorQTL_out/nominal_bpd_risk/`  
Summary Output: `./data/summary/gene_BPD_risk_cell_fraction_interaction.csv` & `./data/summary/gene_BPD_risk_cell_fraction_interaction_summary.csv`  

### MDD Risk SNPs
VCF: `./data/risk_snps/LIBD_maf01_gwas_MDD.vcf.gz`  
Contains 2152/4625 SNPs from:`./data/risk_snps/PGC_depression_genome-wide_significant_makers.txt`  
tensorQTL code: `./code/tensorqtl_risk_interaction_MDD.py`
Raw Output: `./data/tensorQTL_out/nominal_mdd_risk/`  
Summary Output: `./data/summary/gene_MDD_risk_cell_fraction_interaction.csv` & `./data/summary/gene_MDD_risk_cell_fraction_interaction_summary.csv`  


## Genomewide
VCF:  `../genotype_data/mdd_bpd/maf01/mdd_bpd_maf01.vcf.gz` ~11M SNPs

### Nominal - all results
Using `cis.map_nominal()`
tensorQTL code: `./code/tensorqtl_genomewide_nominal_argv.py`  
Raw Output: `./data/tensorQTL_out/genomewide_nominal`  
Summary counts: `./data/summary/cis_genomewide_FDR01.csv`  

### Cis - best SNP per feature
Using `cis.map_cis()`  
tensorQTl code: `./code/tensorqtl_genomewide_cis.py`  
Raw Output: `./data/tensorQTL_out/genomewide_cis`  
Summary counts: TODO  
