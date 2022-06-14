import my_tensorqtl_run
import sys
import os
from tensorqtl import genotypeio, cis, trans
import pandas as pd

out_path = "../data/tensorQTL_out/risk_nominal_mdd/"

for region in ["Amygdala", "sACC"]:
    express, covar = my_tensorqtl_run.get_input_paths("gene", region)
    plink = "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/mdd_bpd/maf01/mdd_bpd_maf01"
    genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df = my_tensorqtl_run.load_data(plink, express, covar, add_chr = True)

    cell_fraction = pd.read_csv("../data/tensorQTL_input/interaction/cell_fraction_"+ region +".csv", index_col = 0)
    
    for cell_type in cell_fraction.columns:
        tag = 'gene_' + region +'_' + cell_type 
        print(tag)
        
        cf_interaction = pd.Series(cell_fraction[cell_type])
        
        nominal_out = cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
        prefix = tag, covariates_df=covariates_df, output_dir = out_path, 
        interaction_s = cf_interaction, maf_threshold_interaction = 0.05, group_s = None, window = 500000, 
        run_eigenmt = False, write_stats = True)
        
    print("DONE " + region)

