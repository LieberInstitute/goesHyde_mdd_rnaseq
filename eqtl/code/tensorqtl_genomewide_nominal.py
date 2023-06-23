import my_tensorqtl_run
from tensorqtl import cis
import sys
import pandas as pd

pair = sys.argv[1].split("_")
feature = pair[0]
region = pair[1]

print("FEATURE = " + feature)
print("REGION = " + region)

expres, covar = my_tensorqtl_run.get_input_paths(feature, region)

prefix = '/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_out/genomewide_nominal_052319'
plink = "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/mdd_bpd/maf01/mdd_bpd_maf01"

genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df = my_tensorqtl_run.load_data(plink, expres, covar, add_chr = True, fix_geno_names = True)

tag = prefix +"/" + feature + "_" + region
print('\n**** STARTING tensorQTL ****')
print("Saving output to: " + tag)

cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix = tag, covariates_df= covariates_df,
                maf_threshold=0.05, interaction_s=None, maf_threshold_interaction=0, 
                group_s=None, window=500000, run_eigenmt=True, output_dir= tag, write_top=False, verbose=False)

print(" DONE !")

