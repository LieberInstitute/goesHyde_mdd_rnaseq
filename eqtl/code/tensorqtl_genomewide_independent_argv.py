import my_tensorqtl_run
from tensorqtl import cis
import pandas as pd
import os.path
import sys

# define paths to data

pair = sys.argv[1].split("_")
region = pair[0]
chrom = pair[1]

print("REGION = " + region)
print("CHR = " + chrom)

expres, covar = my_tensorqtl_run.get_input_paths("gene", region)
prefix = '/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_out/genomewide_independent'
plink = "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/mdd_bpd/maf01/mdd_bpd_maf01"

cis_out_fn = '/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_out/genomewide_cis/gene_' + region + ".csv"
cis_out = pd.read_csv(cis_out_fn, sep=',', index_col=0)

## filter to single chrom
chr_filter = [s.split(':')[0] == chrom for s in cis_out.variant_id]
cis_out = cis_out[chr_filter]
cis_out.shape

genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df = my_tensorqtl_run.load_data(plink, expres, covar, add_chr = True)

tag = prefix +"/ind_gene_" + region + "_" + chrom + "gpu.csv"

## Check dimensions
print('**** STARTING tensorQTL ****')
print("Saving output to: " + tag)

ind_out = cis.map_independent(genotype_df = genotype_df, variant_df = variant_df, cis_df = cis_out,
         phenotype_df = phenotype_df, phenotype_pos_df = phenotype_pos_df, covariates_df = covariates_df,
         group_s=None, maf_threshold=0, fdr=0.01, nperm=10000,window=500000, random_tiebreak=False, 
         logger=None, seed=121,verbose=True)

ind_out.to_csv(tag)

print("**** DONE ****")

