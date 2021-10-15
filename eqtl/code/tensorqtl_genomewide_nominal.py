import my_tensorqtl_run
from tensorqtl import cis
import sys
import pandas as pd
from glob import glob

# define paths to data

for feature in ["gene", "exon", "jxn", "tx"]:
    for region in ["Amygdala", "sACC"]:
        plink, expres, covar = my_tensorqtl_run.get_input_paths(feature, region)
        
        prefix = '/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_out/cis_genomewide_nominal'
        plink = "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/mdd_bpd/maf01/mdd_bpd_maf01"
        
        genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df = my_tensorqtl_run.load_data(plink, expres, covar, add_chr = True)

        tag = prefix +"/" + feature + "_" + region
        print("Saving output to: " + tag)

        cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix = tag, covariates_df= covariates_df,
                           maf_threshold=0, interaction_s=None, maf_threshold_interaction=0, 
                           group_s=None, window=500000, run_eigenmt=True, output_dir= tag, write_top=False, verbose=False)

print(" DONE !")

