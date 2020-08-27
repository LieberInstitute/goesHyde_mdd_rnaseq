data
========

### Summary
* [`merge_drop_swap.R`](merge_drop_swap.R): Combine MDD (n634) and BPD (n540) results in intermediate n1174 samples. Swap BrNums and drop bad rna-samples. Assign genotype match, and drop problem Brains. The result is 1140 sampels this will be the data uploaded to Synpase.
* [`clean_data.R`](clean_data.R): Drop samples based off of csv created in [`ppt_plots.R`](ppt_plots.R) 
