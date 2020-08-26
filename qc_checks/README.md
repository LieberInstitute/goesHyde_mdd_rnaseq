qc_checks
========

### Summary
* [`ppt_plots.R`](ppt_plots.R): Filter samples for totalAssignedGene, overallMapRate, numReads, regional gene expression, and race. Creates plots and qc_dropping_results.csv. 

### [`ppt_plots.R`](ppt_plots.R)

Related files:

Details:

* Filter rna samples for totalAssignedGene, overallMapRate, numReads
* Create sactter plots [`RIN_check_predrop.pdf`](RIN_check_predrop.pdf) & [`RIN_check_postdrop.pdf`](RIN_check_postdrop.pdf)
* Estimate brain regions for each sample, flag mismatches, create boxplot [`region_check_100.pdf`](region_check_100.pdf)
* Flag non-white samples
* Create PCA plots [`pca_log2Rpkm_PC1_2.pdf`](pca_log2Rpkm_PC1_2.pdf) & [`pca_log2Rpkm_PC1_2_combined_datasets.pdf`](pca_log2Rpkm_PC1_2_combined_datasets.pdf)