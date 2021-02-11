#!/bin/sh

mkdir logs/

qsub filter_snps_sacc.sh

qsub filter_snps_amyg.sh
