import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
print(f'PyTorch {torch.__version__}')
print(f'Pandas {pd.__version__}')

def my_tensorqtl_run(plink_prefix_path, expression_bed, covariates_file, prefix):
    print("Starting " + prefix)
    
    # load phenotypes and covariates
    print("Reading Expression files")
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
    covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

    # PLINK reader for genotypes
    print("Reading Plink files")
    pr = genotypeio.PlinkReader(plink_prefix_path)
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
    
    ## Fix genoSample names
    fam = pd.read_table(plink_prefix_path + '.fam', delimiter = " ", names = ["V" + str(i) for i in range(6)])
    genoSamples = [str(v1) + "_" + str(v0) for (v0,v1) in zip(fam.V1, fam.V0)]
    genotype_df.columns = genoSamples
    
    ## Fix chr names (maybe fix in plink?)
    variant_df.chrom = "chr" + variant_df.chrom

    # Filter expression to chrom in snp data
    my_chrom = set(variant_df.chrom)
    chrom_filter = phenotype_pos_df.chr.isin(my_chrom)

    phenotype_df_filter = phenotype_df[chrom_filter]
    phenotype_pos_df_filter = phenotype_pos_df[chrom_filter]
    
    print("Run tensorQTL")
    cis_df = cis.map_cis(
        genotype_df, 
        variant_df, 
        phenotype_df_filter,
        phenotype_pos_df_filter,
        covariates_df=covariates_df,
        seed=9172021)
        
    cis_df.to_csv(prefix + ".csv")

