# PredictDB-LIBD

## 05/03/2021

I'm going to be honest, it's a bit confusing, but I'm going to try to catalog my thoughts and basic understanding of this PredictDB pipeline in here.

There are a few relevant repos:

* https://github.com/hakyimlab/PredictDB_Pipeline_GTEx_v7
* https://github.com/hakyimlab/PredictDB-Tutorial
* https://github.com/hakyimlab/PredictDBPipeline

There is also a write-up in a Google Groups page:

https://groups.google.com/g/predixcanmetaxcan/c/TkBxYkUpNGw/m/Q_mMApRtCQAJ

which is regarding this repo: https://github.com/hakyimlab/PredictDB_Pipeline_GTEx_v7

Hard-coded paths to change:
```
15:55 PredictDB_Pipeline_GTEx_v7 $ grep -R "/group/im-lab/nas40t2/scott/gtex_v7_imputed_europeans/model_training/scripts/" model_training/scripts/
model_training/scripts/whole_blood_test.pbs:#PBS -d /group/im-lab/nas40t2/scott/gtex_v7_imputed_europeans/model_training/scripts/
model_training/scripts/Whole_Blood_test.R:setwd("/group/im-lab/nas40t2/scott/gtex_v7_imputed_europeans/model_training/scripts/")
model_training/scripts/gtex_tiss_chr_elasticnet.pbs:#PBS -d /group/im-lab/nas40t2/scott/gtex_v7_imputed_europeans/model_training/scripts/
model_training/scripts/gtex_tiss_chrom_training.R:setwd("/group/im-lab/nas40t2/scott/gtex_v7_imputed_europeans/model_training/scripts/")
```

There are also a few inputs:
* Expression and covariate data files
* Gene annotation
* Genotype data

Expression and covariate data:
```
NAME              IND1  IND2  IND3 ...
ENSG00000227232.5 -0.81 -0.29 0.31 ...
```

Gene annotation (which might be the `rse_gene` file):
```
chr gene_id           gene_name     start end    gene_type
1   ENSG00000243485.5 MIR1302-2HG   29554 31109  lincRNA
1   ENSG00000237613.2 FAM138A       34554 36081  lincRNA
1   ENSG00000186092.4 OR4F5         69091 70008  protein_coding
1   ENSG00000238009.6 RP11-34P13.7  89295 133723 lincRNA
```

> It is expected at `prepare_data/expression/gencode.v19.genes.patched_contigs.parsed.txt` as we used GENCODE release 19. If you want to use a different one, you should search in the code for the file name or the path stem and change it accordingly.

We used psychENCODE so we should definitely change it.


Genotype data:
```
varID           IND1 IND2  IND3 ...
1_54421_A_G_b37 1    0     0    ...
```

> There must be companion variant annotation files. They are text files with the following format:

```
chromosome pos    varID            ref_vcf alt_vcf R2                 MAF     rsid      rsid_dbSNP150

1          566875 1_566875_C_T_b37 C       T       0.9747600000000001 0.03085 rs2185539 rs2185539
```

## 05/05/2021

Again, for this pipeline, it looks like 
