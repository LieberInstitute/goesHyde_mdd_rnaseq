#!/usr/bin/env python
import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis
import time
import sys
print(f'PyTorch {torch.__version__}')
print(f'Pandas {pd.__version__}')

pair = sys.argv[1].split("_")
feature = pair[0]
region = pair[1]
print("FEATURE = " + feature)
print("REGION = " + region)

outdir='out/genomewide_nominal_062423'
tag = feature + "_" + region + '_gwnom'    # output file prefix

#expres, covar = get_input_paths(feature, region)
expres = 'expression_bed/' + feature + '_' + region + '.bed.gz'
covar = 'covariates_txt/covariates_' + feature + '_' + region + '.txt'

plink = "/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/goesHyde_mdd_rnaseq/genotype_data/mdd_bpd/maf01/mdd_bpd_maf01"

#genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df = load_data(plink, expres, covar, add_chr = True, fix_geno_names = True)
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expres)
print("Phenotype dimensions:", end='')
print(phenotype_df.shape)
covariates_df = pd.read_csv(covar, sep='\t', index_col=0).T
# PLINK reader for genotypes
#print("Reading Plink files: " + plink_prefix_path)
pr = genotypeio.PlinkReader(plink)
#print("Loading Genotypes...", end='')
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
print("Genotype dimensions:", end='')
print(genotype_df.shape)

## Fix genoSample names
#print("Using fam to assign genoSample names")
fam = pd.read_table(plink + '.fam', delimiter = " ", names = ["V" + str(i) for i in range(6)])
genotype_df.columns = [str(v1) + "_" + str(v0) for (v0,v1) in zip(fam.V1, fam.V0)]
# Fix chr names (maybe fix in plink?)
#print("Adding 'chr' to genotype positions")
variant_df.chrom = [s.split(':')[0] for s in list(variant_df.index)]
# Filter expression to chrom in snp data
variant_chrom = set(variant_df.chrom)
express_chrom = set(phenotype_pos_df.chr)
if express_chrom - variant_chrom:
    chrom_filter = phenotype_pos_df.chr.isin(variant_chrom)
    phenotype_df = phenotype_df[chrom_filter]
    phenotype_pos_df = phenotype_pos_df[chrom_filter]

print('\n**** STARTING tensorQTL ****')
print("Saving output to: " + outdir + '/' + tag)

cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix = tag, covariates_df= covariates_df,
                maf_threshold=0.05, window=500000, run_eigenmt=True, output_dir= outdir, write_top=False, verbose=False)

print(" DONE !")

