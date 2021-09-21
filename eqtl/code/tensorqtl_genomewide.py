from my_tensorqtl_run import my_tensorqtl_run
import os.path

# define paths to data

def get_input_paths(feature, region):
    plink_prefix_path = "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/mdd_bpd/maf01/mdd_bpd_maf01"
    expression_bed = '../data/expression_bed/' + feature + '_' + region + '.bed.gz'
    covariates_file = '../data/covariates_txt/covariates_' + feature + '_' + region + '.txt'
    prefix = '/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/eqtl/data/tensorQTL_out/' + feature + '_' + region + '_cis_genomewide'
    return plink_prefix_path, expression_bed, covariates_file, prefix


plink, expres, covar, prefix = get_input_paths("gene", "Amygdala")
my_tensorqtl_run(plink, expres, covar, prefix)



## for feature in ["gene", "exon","jxn","tx"]:
##    for region in ["Amygdala", "sACC"]:
##         plink, expres, covar, prefix = get_input_paths(feature, region)
##         print(expres + " - " + str(os.path.exists(expres)))
##         my_tensorqtl_run(plink, expres, covar, prefix)


