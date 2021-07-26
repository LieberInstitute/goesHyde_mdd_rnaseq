# PredictDB-LIBD

Tutorial link: https://groups.google.com/g/predixcanmetaxcan/c/TkBxYkUpNGw/m/Q_mMApRtCQAJ

Updated repo link: https://github.com/hakyimlab/PredictDB-Tutorial

Inputs for the PredictDB Pipeline:

Gene Annotation file: `/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf`

Genotype: `/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/goesHyde_bipolarMdd_Genotype.rda`

Genotype Covariates: `/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/goesHyde_bipolarMdd_Genotypes_mds.csv`

TMM-Normalized Expression: `/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/goesHyde_mdd_rnaseq_Amygdala.expression.bed.gz`

Combined Covariates: `/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/goesHyde_mdd_rnaseq_Amygdala.combined_covariates.txt`

```
There must be companion variant annotation files. They are text files with the following format:

chromosome pos    varID            ref_vcf alt_vcf R2                 MAF     rsid      rsid_dbSNP150

1          566875 1_566875_C_T_b37 C       T       0.9747600000000001 0.03085 rs2185539 rs2185539


The variant ids need not be in the "{chr}_{pos}_{NEA}_{EA}_b37" format; the variant annotation file must provide the mapping from variant id to rsid.


Genotype and annotation files are expected to be split by chromosome.
```
