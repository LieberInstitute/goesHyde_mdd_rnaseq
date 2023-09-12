import my_tensorqtl_run
import sys
import os
from tensorqtl import genotypeio, cis, trans
import pandas as pd
import session_info

out_path = "../data/tensorQTL_out/interaction_nominal/"

for feature in ["tx"]:
#for feature in ["gene", "tx"]:
    for region in ["Amygdala", "sACC"]:
        
        express, covar = my_tensorqtl_run.get_input_paths(feature, region)
        plink = "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/mdd_bpd/maf01/mdd_bpd_maf01"
        genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df = my_tensorqtl_run.load_data(plink, express, covar, add_chr = True, fix_geno_names = True)
        
        cell_fraction = pd.read_csv("../data/tensorQTL_input/interaction/cell_fraction_"+ region +".csv", index_col = 0)
        
        for cell_type in cell_fraction.columns:
            tag = feature '_' + region +'_' + cell_type 
            print(tag)
            
            cf_interaction = cell_fraction[[cell_type]]
            
            nominal_out = cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
            prefix = tag, covariates_df=covariates_df, output_dir = out_path, 
            interaction_df = cf_interaction, maf_threshold_interaction = 0.05, group_s = None, window = 500000, 
            run_eigenmt = False, write_stats = True)
        
        print("DONE " + region)

session_info.show()

