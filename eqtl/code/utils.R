
get_resid_expres <- function(rse, mod, tpm = FALSE){
  mod <- mod[colnames(rse),]
  if(tpm) rpkm <- assays(rse)$tpm else rpkm <- recount::getRPKM(rse, "Length")
  ## residualize expression
  exprs <- log2(rpkm + 1)
  exprs <- cleaningY(exprs, mod, P = 1)
  colnames(exprs) <- rse$genoSample
  return(exprs)
}

