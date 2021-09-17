import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
print(f'PyTorch {torch.__version__}')
print(f'Pandas {pd.__version__}')

# define paths to data
plink_prefix_path = "LIBD_maf01_gwas_BPD"
expression_bed = 'pheno_data/gene_Amygdala.bed.gz'
covariates_file = 'pheno_data/covariates_gene_Amygdala.txt'
prefix = 'gene_amygdala'

# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

cis.map_nominal(genotype_df, variant_df,
                phenotype_df.loc[phenotype_pos_df['chr']=='chr18'],
                phenotype_pos_df.loc[phenotype_pos_df['chr']=='chr18'],
                prefix, covariates_df=covariates_df)