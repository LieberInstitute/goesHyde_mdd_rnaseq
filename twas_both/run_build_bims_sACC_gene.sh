#!/bin/bash

echo "**** Job starts ****"
date

## Clean up older files
# rm -fr trash

## Could be fancier and use the date or something
mkdir -p trash
mv logs/build_bims_*.txt trash/
rm -r trash/sACC_gene
mv sACC_gene trash/

## Create logs dir if needed
mkdir -p logs
## Submit new job
qsub build_bims_sACC_gene.sh

echo "**** Job ends ****"
date
