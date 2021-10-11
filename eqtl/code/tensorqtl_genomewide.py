import my_tensorqtl_run
from tensorqtl import cis
import pandas as pd
import os.path

# define paths to data

for region in ["Amygdala", "sACC"]:
    plink, expres, covar, prefix = my_tensorqtl_run.get_input_paths("gene", region)
    print(prefix)
    plink = "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/mdd_bpd/maf01/mdd_bpd_maf01"
    genotype_df, variant_df, phenotype_df_filter, phenotype_pos_df_filter, covariates_df = my_tensorqtl_run.load_data(plink, expres, covar)
    
    cis_out = cis.map_cis(genotype_df, variant_df, phenotype_df_filter, phenotype_pos_df_filter, covariates_df=covariates_df,
                group_s=None, maf_threshold=0, beta_approx=True, nperm=10000,
                window=500000, random_tiebreak=False, logger=None, seed=None,
               verbose=True, warn_monomorphic=False)

    cis_out.to_csv(prefix + ".csv")

