
get_resid_expres <- function(rse, mod, tpm = FALSE) {
    mod <- mod[colnames(rse), ]
    if (tpm) rpkm <- assays(rse)$tpm else rpkm <- recount::getRPKM(rse, "Length")
    ## residualize expression
    exprs <- log2(rpkm + 1)
    exprs <- cleaningY(exprs, mod, P = 1)
    colnames(exprs) <- rse$genoSample
    return(exprs)
}

plot_my_eqtl <- function(genotype, expression, annotation){
    return(NULL)
}


read_adj_filter <- function(parquet_files, cutoff = 0.05, pval_name = "pval_nominal") {
    eqtl_out <- do.call("rbind", map(parquet_files, parquet_read)) %>%
        mutate(FDR = p.adjust(get(pval_name), "fdr"))
    message("n pairs: ", nrow(eqtl_out))
    # filter
    eqtl_out <- eqtl_out %>%
        filter(FDR < cutoff)
    message("n pairs FDR<", cutoff, ": ", nrow(eqtl_out))
    return(eqtl_out)
}


