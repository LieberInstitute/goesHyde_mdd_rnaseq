
rse_to_bed <- function(rse) {
    rr_df <- as.data.frame(rowRanges(rse)) %>%
        tibble::rownames_to_column("ID") %>%
        dplyr::arrange(seqnames, start) %>%
        dplyr::mutate(end = start + 1) %>%
        dplyr::select(`#Chr` = seqnames, start, end, ID)


    counts <- SummarizedExperiment::assays(rse)$counts
    colnames(counts) <- rse$genoSample
    # counts = round(counts, 3) ## round to 3 ?

    counts <- as.data.frame(counts) %>%
        tibble::rownames_to_column("ID") %>%
        dplyr::left_join(rr_df, ., by = "ID")

    return(counts)
}
