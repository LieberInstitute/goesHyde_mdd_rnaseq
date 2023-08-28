#' Convert assay from a SummarizedExperiment to tensorQTL compatible bed file
#'
#' @param rse 
#' @param assay_name
#' @return A `data.frame()` that can be saved with `write.table(sep = "\t",quote = FALSE, row.names = FALSE)` to import to tensorQTL
#' @export
#'
#' @examples bed <- rse_to_bed(rse_gene)
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr arrange mutate select left_join
#' @importFrom SummarizedExperiment assays
rse_to_bed <- function(rse, assay_name = "logcounts") {
  
  stopifnot(assay_name %in% assayNames(rse))
  
  rr_df <- as.data.frame(rowRanges(rse)) |>
    tibble::rownames_to_column("ID") |>
    dplyr::mutate(start = ifelse(strand == "+", start, end),
                  end = start + 1) |>
    dplyr::arrange(seqnames, start) |>
    dplyr::select(`#Chr` = seqnames, start, end, ID)
  
  message(Sys.time(), " - Converting assay '", assay_name, "' to bed format: (", nrow(rse), ", ", ncol(rse),")")
  counts <- SummarizedExperiment::assays(rse)[[assay_name]]
  colnames(counts) <- rse$genoSample
  
  counts <- rr_df |> 
    dplyr::left_join(as.data.frame(counts) |> 
                       tibble::rownames_to_column("ID"), 
                     by = "ID")
  return(counts)
}
