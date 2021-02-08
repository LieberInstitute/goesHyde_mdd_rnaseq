
run_DE <- function(rse, model, run_voom = TRUE){
  
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
    fit = lmFit(log_tm, model)
  }
  
  ## limma
  
  eBayes_out = eBayes(fit)
  
  #### because of more than 1 component, computing F statistics instead of t=-statistics
  topTable_out = topTable(eBayes_out, coef=2:3,
                          p.value = 1, number=nrow(rse))
  #### reorderoing based on genes in rse_gene
  topTable_out = topTable_out[rownames(rse),]
  
  ## significance levels EXTRACT INDIVIDUAL COMPARISON P-VALUES THAT ARE NOT IN TOP TABLE
  pvalMat = as.matrix(eBayes_out$p.value)[,2:3]
  
  ### check top p-values
  head(pvalMat[order(pvalMat[,"PrimaryDxControl"]), ])
  
  qvalMat = pvalMat
  qvalMat[,1:2] = p.adjust(pvalMat[,1:2],method="fdr")
  colnames(pvalMat) = paste0("P_",colnames(pvalMat))
  colnames(qvalMat) = paste0("q_",colnames(qvalMat))
  
  topTable_out = cbind(topTable_out,cbind(pvalMat, qvalMat))
  head(topTable_out)
  sum(topTable_out$q_PrimaryDxControl < 0.05)
  
  sig05 <- sum(topTable_out$q_PrimaryDxControl < 0.05)
  sig01 <- sum(topTable_out$q_PrimaryDxControl < 0.01)

  message(paste("PrimaryDxControl < 0.05:", sig05))
  message(paste("PrimaryDxControl < 0.01:", sig01))
  message("PrimaryDxBipolar < 0.05: ",
          sum(topTable_out$q_PrimaryDxBipolar < 0.05))

  message("PrimaryDxBipolar < 0.01: ",
          sum(topTable_out$q_PrimaryDxBipolar < 0.01))

  return(topTable_out)
}
