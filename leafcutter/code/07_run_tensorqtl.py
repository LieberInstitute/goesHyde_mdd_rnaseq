import my_tensorqtl_run
from tensorqtl import cis
import sys
import pandas as pd


region = sys.argv[1]

print("Loading data for: " + region)

prefix = '/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/leafcutter/data/tensorQTL_out'
plink = "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/mdd_bpd/maf01/mdd_bpd_maf01"


## splice.bed.gz
expres = "../data/qqnorm/qqnorm_" + region + ".bed.gz"
covar = "../data/covariates/covariates_" + region + ".txt"

genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df = my_tensorqtl_run.load_data(plink, expres, covar, add_chr = True, fix_geno_names = True)


tag = prefix +"/LC-clusters_" + region
print('\n**** STARTING tensorQTL ****')
print("Saving output to: " + tag)

cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix = tag, covariates_df= covariates_df,
                maf_threshold=0.05, interaction_s=None, maf_threshold_interaction=0, 
                group_s=None, window=500000, run_eigenmt=True, output_dir= tag, write_top=False, verbose=False)

print(" DONE !")