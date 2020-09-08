#!/bin/bash

module load plink
plink --vcf topmed_mdd_602sample_090120_maf005.vcf.gz --make-bed --out topmed_mdd_602sample_090120_maf005
