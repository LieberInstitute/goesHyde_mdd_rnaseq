# PredictDB-LIBD

# PredictDB-LIBD Workflow


TODO make a DAG representing this pipeline but with more detail, include all inputs/outputs.

## Get Normalized Expression
The names of these scripts are a bit misleading because it says "inverse quantile normalized" but really it's, by default, trimmed mean of M (TMM) normalized. One of the authors of the PredictDB pipeline has said that the normalization method in particular does not matter all too much, so I have kept it TMM.

Scripts:
- `predixcan_pipeline/code/01_get_inv_quantile_norm/01_get_inv_quantile_norm.sh`
  - First script of PredictDB pipeline
  - `qsub`
- `predixcan_pipeline/code/01_get_inv_quantile_norm/01_prepare_gene_expression.R`
  - Inputs:
    - `exprs_cutoff/rse_gene.Rdata`
      - Expression file
    - `genotype_data/topmed_mdd_602sample_090120_maf005.bim/bed/fam`
      - Input genotype
  - Outputs:
    - `predixcan_pipeline/processed-data/01_get_inv_quantile_norm/{region}_rse_sub_sample_IDs.txt`
      - Sample IDs of RNA-Seq from one brain region
    - `predixcan_pipeline/processed-data/01_get_inv_quantile_norm/topmed_mdd_602sample_090120_maf005_Amygdala_sorted.bim/bed/fam`
      - Genotypes of RNA-Seq samples above
      - From one brain region
    - `predixcan_pipeline/processed-data/01_get_inv_quantile_norm/topmed_mdd_602sample_090120_maf005_Amygdala_sorted.vcf.gz`
      - Same as above in VCF format
    - `predixcan_pipeline/processed-data/01_get_inv_quantile_norm/{region}_vcf_samples.txt`
      - Sample IDs extracted again
    - `predixcan_pipeline/processed-data/01_get_inv_quantile_norm/{region}_gene_counts.gct`
      - Gene counts file necessary for `eqtl_prepare_expression.py`
    - `predixcan_pipeline/processed-data/01_get_inv_quantile_norm/{region}_tpm.gct`
      - Gene counts file necessary for `eqtl_prepare_expression.py`
    - `predixcan_pipeline/processed-data/01_get_inv_quantile_norm/{region}_samp_part_lookup.txt`
      - Sample/participant lookup table
      - Two identical columns
- `predixcan_pipeline/code/01_get_inv_quantile_norm/KJ_eqtl_prepare_expression.py`
  - KJ Benjamin's custom version of GTEx's `eqtl_prepare_expression.py` script
  - Inputs:
    - `predixcan_pipeline/processed-data/01_get_inv_quantile_norm/{region}_tpm.gct`
    - `predixcan_pipeline/processed-data/01_get_inv_quantile_norm/{region}_gene_counts.gct`
    - `/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf`
      - Gene annotation
    - `predixcan_pipeline/processed-data/01_get_inv_quantile_norm/{region}_samp_part_lookup.txt`
    - `/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/vcf_chr_list.txt`
      - List of chromosomes
    - `prefix_Amygdala`
      - Project naming scheme
      - `goesHyde_mdd_rnaseq_Amygdala`
    - Outputs:
      - `${prefix}.expression.bed.gz`
      - `${prefix}.expression.bed.gz.tbi`

### Get PEER Covariates

We next use a [conda environment called `r-peer`](https://anaconda.org/bioconda/r-peer) to call PEER and calculate covariates. `run_PEER.R` is from the [broadinstitute/gtexpipeline](https://github.com/broadinstitute/gtex-pipeline/blob/75feccdf238f6fc69c3a8e5aed57bae7fb3be953/qtl/src/run_PEER.R) repository, reproduced here for ease of access.

Scripts:
- `predixcan_pipeline/code/01_get_inv_quantile_norm/02_run_PEER.sh`
  - `qsub`
- `predixcan_pipeline/code/01_get_inv_quantile_norm/run_PEER.R`
  - Inputs:
    - `num_peer`
      - Variable representing number of PEER factors [based on recommendation](https://github.com/broadinstitute/gtex-pipeline/tree/75feccdf238f6fc69c3a8e5aed57bae7fb3be953/qtl#2-calculate-peer-factors), `60`
    - `prefix_Amygdala`
      - Project naming scheme
      - `goesHyde_mdd_rnaseq_Amygdala`
    - `predixcan_pipeline/code/01_get_inv_quantile_norm/run_PEER.R processed-data/01_get_inv_quantile_norm/${prefix_Amygdala}.expression.bed.gz`
      - Expression VCF generated in previous step
  - Outputs:
    - `${prefix_Amygdala}.PEER_residuals.txt`
    - `${prefix_Amygdala}.PEER_alpha.txt`
    - `${prefix_Amygdala}.PEER_covariates.txt`

### Combine All Covariates

Simple script designed to combine all relevant into one file.

Scripts:
- `predixcan_pipeline/code/01_get_inv_quantile_norm/03_combine_covs.sh`
  - `qsub`
- `predixcan_pipeline/code/01_get_inv_quantile_norm/03_convert_rda.R`
  - Converts the covariates in the genotype file into a text file
  - Reformats the covariates output by PEER as csv
  - Inputs:
    - `/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/goesHyde_bipolarMdd_Genotypes_mds.rda`
    - `exprs_cutoff/rse_gene.Rdata`
    - `predixcan_pipeline/processed-data/01_get_inv_quantile_norm/goesHyde_mdd_rnaseq_Amygdala.PEER_covariates.txt`
  - Outputs:
    - `predixcan_pipeline/processed-data/01_get_inv_quantile_norm/goesHyde_bipolarMdd_Genotypes_mds.csv`
    - `predixcan_pipeline/processed-data/01_get_inv_quantile_norm/goesHyde_mdd_rnaseq_Amygdala.PEER_covariates_REFORMATTED.txt`
      - Doesn't seem like I actually used this output
- `predixcan_pipeline/code/01_get_inv_quantile_norm/combine_covariates.py`
  - Inputs:
    - `prefix_Amygdala`
      - `goesHyde_mdd_rnaseq_Amygdala`
    - `genotype_pcs`
      - Extracted from `/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/genotype_data/goesHyde_bipolarMdd_Genotypes_mds.rda`
      - `/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/goesHyde_bipolarMdd_Genotypes_mds.csv`
    - `processed-data/01_get_inv_quantile_norm/${prefix_Amygdala}.PEER_covariates.txt`
      - PEER covariates from previous step

## Prep Inputs

### Custom Genotype

As I may have mentioned earlier, I ran into a problem in the course of constructing this pipeline where the genotype file does not contain the full IDs of the samples, and only represents the "family" ID. In order to remedy this, Josh Stolz directed me towards the script responsible for generating the genotype file and I removed the lines which truncated the full sample ID.

- Scripts:
  - `predixcan_pipeline/code/02_prep_inputs/00_call_pull_genotype_data.sh`
  - `predixcan_pipeline/code/02_prep_inputs/01_pull_genotype_data.R`
    - Modified script to contain the full sample IDs of all of the samples
    - Does not output `mds` file like original scirpt
    - Inputs:
      - `predixcan_pipeline/data/rse_gene_GoesZandi.rda`
      - `genotype_data/goesHyde_mdd_Genotypes_maf01_geno10_hwe1e6.bim/bed/fam`
      - `/dcs01/ajaffe/Annotation/dbsnp142_common.txt`
      - `genotype_data/dbSNP.Rdata`
      - `/dcl01/lieber/ajaffe/Brain/Imputation/Merged/hg19ToHg38.over.chain`
    - Outputs:
      - `predixcan_pipeline/processed-data/02_prep_inputs/goesHyde_bipolarMdd_Genotypes_PredictDB_NO-MDS.rda`


### Prep SNP Annotation and Genotype

  - `predixcan_pipeline/code/02_prep_inputs/02_prep_predictdb.sh`
    - Calls 4 scripts, 2 written by me and 2 in `PredictDB-Tutorial`
  - `predixcan_pipeline/PredictDB-Tutorial/code/parse_gtf.py`
    - Inputs:
      - `/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/GENCODE/GRCh38_hg38/gencode.v25.annotationGRCh38.gtf`
    - Outputs:
      - `/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/02_prep_inputs/gencode.v25.annotationGRCh38_PARSED.gtf`
  - `predixcan_pipeline/code/02_prep_inputs/03_process_snp_anno.R`
    - Takes the genotype and SNP annotation and formats them for compatibility with the `PredictDB-Tutorial` elastic net script
    - Removes SNPs with `NA` RSIDs and creates ad hoc genomic variants in hg38 for cross-reference purposes only
    - Saves genotype file by chromosome
    - Inputs:
      - `predixcan_pipeline/processed-data/02_prep_inputs/goesHyde_bipolarMdd_Genotypes_PredictDB_NO-MDS.rda`
    - Outputs:
      - `predixcan_pipeline/processed-data/02_prep_inputs/split_geno/split_snp_geno.chr{1-22}.txt`
      - `predixcan_pipeline/processed-data/02_prep_inputs/snp_annot_prep.txt`

### Prep Normalized Expression File

  - `predixcan_pipeline/code/02_prep_inputs/04_process_norm_exp.R`
    - Essentially just transposes the expression table for compatibility
    - Inputs:
      - `predixcan_pipeline/processed-data/01_get_inv_quantile_norm/goesHyde_mdd_rnaseq_Amygdala.expression.bed.gz`
    - Outputs:
      - `predixcan_pipeline/processed-data/02_prep_inputs/transformed_expression.txt`

## Compute Weights

Last of the preliminary steps. Here I simply call their elastic net model handler.

Scripts:
  - `predixcan_pipeline/code/03_compute_weights/01_compute_weights.sh`
  - `predixcan_pipeline/PredictDB-Tutorial/code/gtex_tiss_chrom_training.R`
    - Handles inputs and calls `main` for the elastic net model
      - `predixcan_pipeline/PredictDB-Tutorial/code/gtex_v7_nested_cv_elnet.R`
    - Inputs:
        - Chromosome `1-22`
        - `snp_annot_file <- "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/02_prep_inputs/snp_annot/snp_annot.chr" %&% chrom %&% ".txt"`
        - `gene_annot_file <- "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/02_prep_inputs/gencode.v25.annotationGRCh38_PARSED.gtf"`
        - `genotype_file <- "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/02_prep_inputs/split_geno/split_snp_geno.chr" %&% chrom %&% ".txt"`
        - `expression_file <- "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/02_prep_inputs/transformed_expression.txt"`
        - `covariates_file <- "/dcl01/lieber/ajaffe/lab/goesHyde_mdd_rnaseq/predixcan_pipeline/processed-data/01_get_inv_quantile_norm/goesHyde_mdd_rnaseq_Amygdala.combined_covariates.txt"`
        - `prefix <- "goesHyde_mdd_Amygada_Model_training"`
      - Outputs:
        - `predixcan_pipeline/summary/goesHyde_mdd_Amygada_Model_training_chr{1-22}_model_summaries.txt`
        - `predixcan_pipeline/summary/goesHyde_mdd_Amygada_Model_training_chr{1-22}_summary.txt`
        - `predixcan_pipeline/weights/goesHyde_mdd_Amygada_Model_training_chr{1-22}_weights.txt`
