#!/bin/bash

echo "**** Job starts ****"
date

## Clean up older files
# rm -fr trash

## Could be fancier and use the date or something
mkdir -p trash
mv logs/*.txt trash/
rm -r trash/processed-data
mv Amygdala_gene processed-data/

## Create logs dir if needed
mkdir -p logs
mkdir -p processed-data/01_get_inv_quantile_norm
mkdir -p processed-data/02_prep_inputs  

echo "**** Job ends ****"
date
