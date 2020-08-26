data
========

### Summary
* [`merge_gene_raw_BPD_MDD.R`](merge_gene_raw_BPD_MDD.R): Combine MDD (n634) and BPD (n540), no filtering, results in n1174 sample object
* [`mdd_swap.R`](mdd_swap.R): Swap BrNums and drop bad rna-samples. Assign genotype match, and drop problem Brains. The result is 1140 sampels this will be the data uploaded to Synpase.
* [`clean_data.R`](clean_data.R): Drop samples based off of csv created in [`ppt_plots.R`](ppt_plots.R) 
