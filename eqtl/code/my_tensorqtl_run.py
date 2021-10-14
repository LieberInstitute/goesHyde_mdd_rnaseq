import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
print(f'PyTorch {torch.__version__}')
print(f'Pandas {pd.__version__}')

def load_data(plink_prefix_path, expression_bed, covariates_file, add_chr = False):
    
    # load phenotypes and covariates
    print("Reading Expression files: " + expression_bed)
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
    print("Phenotype dimensions:")
    print(phenotype_df.shape)
    
    covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T
    
    print("." *20 )
    # PLINK reader for genotypes
    print("Reading Plink files: " + plink_prefix_path)
    pr = genotypeio.PlinkReader(plink_prefix_path)
    print("Loading Genotypes...", end='')
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
    print("Genotype dimensions:", end='')
    print(genotype_df.shape)
    
    ## Fix genoSample names
    fam = pd.read_table(plink_prefix_path + '.fam', delimiter = " ", names = ["V" + str(i) for i in range(6)])
    genoSamples = [str(v1) + "_" + str(v0) for (v0,v1) in zip(fam.V1, fam.V0)]
    genotype_df.columns = genoSamples
    
    ## Fix chr names (maybe fix in plink?)
    if add_chr:
        print("Adding 'chr' to genotype positions")
        variant_df.chrom = "chr" + variant_df.chrom
        variant_df.replace("chr23", "chrX")

    # Filter expression to chrom in snp data
    my_chrom = set(variant_df.chrom)
    chrom_filter = phenotype_pos_df.chr.isin(my_chrom)
    print("Chromosomes present in geno: " + str(len(my_chrom)))

    phenotype_df_filter = phenotype_df[chrom_filter]
    phenotype_pos_df_filter = phenotype_pos_df[chrom_filter]
    
    print("filter to " + str(phenotype_df_filter.shape[0]) + "/" + str(phenotype_df.shape[0]) + " features")
    
    return(genotype_df, variant_df, phenotype_df_filter, phenotype_pos_df_filter, covariates_df)
    
def get_input_paths(feature, region):
    plink_prefix_path = "../data/risk_snps/LIBD_maf01_gwas_BPD_" + region
    expression_bed = '../data/expression_bed/' + feature + '_' + region + '.bed.gz'
    covariates_file = '../data/covariates_txt/covariates_' + feature + '_' + region + '.txt'
    return plink_prefix_path, expression_bed, covariates_file



