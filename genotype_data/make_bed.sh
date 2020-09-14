#!/bin/bash

module load plink
#plink --vcf topmed_mdd_602sample_090120_maf005.vcf.gz --make-bed --out topmed_mdd_602sample_090120_maf005
plink --vcf Illumina_HumanHap650Yv3_A_mdd.vcf.gz --make-bed --out Illumina_HumanHap650Yv3_A_mdd
plink --vcf Illumina_Omni2.5-8_v1.1_mdd.vcf.gz --make-bed --out Illumina_Omni2.5-8_v1.1_mdd
plink --vcf Illumina_Omni2.5-8_v1.3._mdd.vcf.gz --make-bed --out Illumina_Omni2.5-8_v1.3._mdd
plink --vcf Illumina_Omni2.5-8_v1.5_56_mdd.vcf.gz --make-bed --out Illumina_Omni2.5-8_v1.5_56_mdd
plink --vcf Illumina_Omni2.5_v1.3_8_mdd.vcf.gz --make-bed --out Illumina_Omni2.5_v1.3_8_mdd
plink --vcf Illumina_Omni2.5_v1.4_136_mdd.vcf.gz --make-bed --out Illumina_Omni2.5_v1.4_136_mdd
plink --vcf Illumina_Omni5-Quad_mdd.vcf.gz --make-bed --out Illumina_Omni5-Quad_mdd
