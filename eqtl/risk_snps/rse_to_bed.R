
rse_to_bed <- function(rse){
  rr_df <- as.data.frame(rowRanges(rse)) %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::select(`#chr` = seqnames, left = start, right = end, phenotype_id = gencodeID, group_id = gene_type, strand)
  
    counts <- SummarizedExperiment::assays(rse)$counts
    colnames(counts) <- rse$genoSample
    # counts = round(counts, 3) ## round to 3 ?
    
    counts <- as.data.frame(counts) %>%
      tibble::rownames_to_column("phenotype_id") %>%
      dplyr::left_join(rr_df, ., by = "phenotype_id")
  
    return(counts)
}
