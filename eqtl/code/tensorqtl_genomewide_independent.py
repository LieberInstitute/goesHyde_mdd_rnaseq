import my_tensorqtl_run
from tensorqtl import cis
import pandas as pd
import os.path

# define paths to data

for region in ["Amygdala", "sACC"]:
    expres, covar = my_tensorqtl_run.get_input_paths("gene", region)
    prefix = '/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_out/genomewide_independent'
    plink = "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/mdd_bpd/maf01/mdd_bpd_maf01"
    
    cis_out_fn = '/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_out/genomewide_cis/gene_' + region + ".csv"
    cis_out = pd.read_csv(cis_out_fn, sep=',', index_col=0)
    
    genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df = my_tensorqtl_run.load_data(plink, expres, covar, add_chr = True)
    
    output = prefix +"/gene_" + region + ".csv"
    
    ## Check dimensions
    print("Variant.df shape: ")
    variant_df.shape
    genotype_df.shape[0] == phenotype_df.shape[1]
    
    print('**** STARTING tensorQTL ****')
    print("Saving output to: " + output)
    
    ind_out = cis.map_independent(genotype_df = genotype_df, variant_df = variant_df, cis_df = cis_out,
                phenotype_df = phenotype_df, phenotype_pos_df = phenotype_pos_df, covariates_df = covariates_df,
                group_s=None, maf_threshold=0, nperm=10000,
                window=500000, random_tiebreak=False, logger=None, seed=119,
               verbose=True)

    ind_out.to_csv(output)

