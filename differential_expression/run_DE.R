
run_DE <- function(rse, model, run_voom = TRUE, save_eBayes = FALSE, coef = c("PrimaryDxControl","PrimaryDxBipolar")){
  
  ## limma
  eBayes_out = get_eBayes(rse = rse, model = model, run_voom = run_voom)
  
  #### because of more than 1 component, computing F statistics instead of t=-statistics
  # message("Calc top tables - coefficents:", paste(coef, collpase = " "))
  message("coef in eBayes: ", all(coef %in% colnames(eBayes_out$coefficients)))
  
  topTable_out = topTable(eBayes_out, coef=coef,
                          p.value = 1, number=nrow(rse))
  #### reorderoing based on genes in rse_gene
  topTable_out = topTable_out[rownames(rse),]
  
  ## significance levels EXTRACT INDIVIDUAL COMPARISON P-VALUES THAT ARE NOT IN TOP TABLE
  pvalMat = as.matrix(eBayes_out$p.value)[,coef]
  
  ### check top p-values
  head(pvalMat[order(pvalMat[,"PrimaryDxControl"]), ])
  
  qvalMat = pvalMat
  qvalMat[,1:2] = p.adjust(pvalMat[,1:2],method="fdr")
  colnames(pvalMat) = paste0("P_",colnames(pvalMat))
  colnames(qvalMat) = paste0("q_",colnames(qvalMat))
  
  topTable_out = cbind(topTable_out,cbind(pvalMat, qvalMat))
  
  n_pVal <- report_top_pVal(topTable_out = topTable_out, cols = paste0("q_", coef))
  print(n_pVal)

  if(save_eBayes){
    output <- list(topTable = topTable_out,
                   eBayes = eBayes_out)
    return(output)
  } else {
    return(topTable_out)
  }
  
}

get_eBayes <- function(rse, model, run_voom = TRUE){
  
  if(run_voom){
    message("Calc Norm Factors: ", Sys.time())
    dge_out = DGEList(counts = assays(rse)$counts,
                      genes = rowData(rse))
    dge_norm = calcNormFactors(dge_out)
    voom_out = voom(dge_norm, model, plot=FALSE)
    
    message("Limma: ", Sys.time())
    fit = lmFit(voom_out)
  } else {
    log_tmp = log2(assays(rse)$tpm + 1)
    
    message("Limma: ", Sys.time())
    fit = lmFit(log_tmp, model)
  }
  
  ## limma
  message("eBayes: ", Sys.time())
  eBayes_out = eBayes(fit)
  return(eBayes_out)
}

report_top_pVal <- function(topTable_out, cols){
  
  tt_values <- topTable_out[,cols]
  cutoffs <- list(0.05, 0.01)
  names(cutoffs) <- paste("< ", cutoffs)
  n_sig <- map(cutoffs, ~colSums(tt_values < .x))
  sig_table <- t(do.call("rbind", n_sig))
  return(sig_table)

}
