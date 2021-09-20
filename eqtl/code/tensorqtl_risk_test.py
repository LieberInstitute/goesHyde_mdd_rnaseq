import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
print(f'PyTorch {torch.__version__}')
print(f'Pandas {pd.__version__}')

# define paths to data
plink_prefix_path = "../data/risk_snps/LIBD_maf01_gwas_BPD_Amygdala"
expression_bed = '../data/expression_bed/gene_Amygdala.bed.gz'
covariates_file = '../data/covariates_txt/covariates_gene_Amygdala.txt'
prefix = '/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_out/Amyg_test_cis'

# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

covariates_df.shape
genotype_df.shape
phenotype_df.shape

## Current fixes
genotype_df.columns = phenotype_df.columns
variant_df.chrom = "chr" + variant_df.chrom
variant_df.chrom.value_counts()

my_chrom = set(variant_df.chrom)
chrom_filter = phenotype_pos_df.chr.isin(my_chrom)
phenotype_pos_df.chr in my_chrom

phenotype_df_filter = phenotype_df[chrom_filter]
phenotype_pos_df_filter = phenotype_pos_df[chrom_filter]

cis_df = cis.map_cis(genotype_df, variant_df, 
                     phenotype_df_filter,
                     phenotype_pos_df_filter,
                     covariates_df=covariates_df,
                     seed=9172021)

cis_df.to_csv(prefix + ".csv")

