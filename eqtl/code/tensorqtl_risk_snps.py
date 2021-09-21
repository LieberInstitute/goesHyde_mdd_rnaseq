from my_tensorqtl_run import my_tensorqtl_run
import os.path

# define paths to data

def get_input_paths(feature, region):
    plink_prefix_path = "../data/risk_snps/LIBD_maf01_gwas_BPD_" + region
    expression_bed = '../data/expression_bed/' + feature + '_' + region + '.bed.gz'
    covariates_file = '../data/covariates_txt/covariates_' + feature + '_' + region + '.txt'
    prefix = '/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_out/' + feature + '_' + region + '_cis_risk_bpd'
    return plink_prefix_path, expression_bed, covariates_file, prefix




for feature in ["gene", "exon","jxn","tx"]:
    for region in ["Amygdala", "sACC"]:
         plink, expres, covar, prefix = get_input_paths(feature, region)
         my_tensorqtl_run(plink, expres, covar, prefix)


