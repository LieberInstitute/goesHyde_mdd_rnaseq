import my_tensorqtl_run
from tensorqtl import cis
import pandas as pd
import os.path

# define paths to data

plink, expres, covar, prefix = my_tensorqtl_run.get_input_paths("gene", "Amygdala")
print(prefix)
prefix = prefix + "_nominal"
plink = "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/mdd_bpd/maf01/mdd_bpd_maf01"

genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df = my_tensorqtl_run.load_data(plink, expres, covar, add_chr = True)
    
## cis_out = cis.map_cis(genotype_df, variant_df, phenotype_df_filter, phenotype_pos_df_filter, covariates_df=covariates_df,
## group_s=None, maf_threshold=0, beta_approx=True, nperm=10000,
##                 window=500000, random_tiebreak=False, logger=None, seed=None,
##                verbose=True)

tensor_out = cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix, covariates_df= covariates_df,
maf_threshold=0, interaction_s=None, maf_threshold_interaction=0, group_s=None, window=500000, run_eigenmt=False,
output_dir= prefix +"gene_Amygdala", write_top=False, write_stats=True, logger=None, verbose=True)

tensor_out.to_csv(prefix + ".csv")

