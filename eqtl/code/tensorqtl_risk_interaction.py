import my_tensorqtl_run
import sys
import os
from tensorqtl import genotypeio, cis, trans
import pandas as pd

out_path = "../data/tensorQTL_out/nominal_bpd_risk/"

for region in ["Amygdala", "sACC"]:
    plink, express, covar, prefix = my_tensorqtl_run.get_input_paths("gene", region)
    genotype_df, variant_df, phenotype_df_filter, phenotype_pos_df_filter, covariates_df = my_tensorqtl_run.load_data(plink, express, covar)
    
    cell_fraction = pd.read_csv("../data/interaction/cell_fraction_"+ region +".csv", index_col = 0)
    
    for cell_type in cell_fraction.columns:
        tag = 'gene_' + region +'_' + cell_type 
        print(tag)
        ##out_path_ct = out_path + tag + "/"
        cf_interaction = pd.Series(cell_fraction[cell_type])
        
        nominal_test = cis.map_nominal(genotype_df, variant_df, phenotype_df_filter, phenotype_pos_df_filter,
        prefix = tag, covariates_df=covariates_df, output_dir= out_path + "parquet_data/", 
        interaction_s = cf_interaction, maf_threshold_interaction=0, group_s=None, window=1000000, 
        run_eigenmt=True, write_top = False)
        
        nominal_test.to_csv('../data/tensorQTL_out/nominal_bpd_risk/' + tag + ".csv")

