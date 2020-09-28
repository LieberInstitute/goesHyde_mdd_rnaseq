library(ggplot2)
library(ggbeeswarm)
library(here)

source(here("main_colors.R"))

make_expDx_table <- function(exprs, i, pd){
  exprs_i <- as.data.frame(exprs[i,])
  colnames(exprs_i) <- c("exp")
  exprs_i$Dx <- pd[["Dx"]]
  return(exprs_i)
}

exp_title <- function(out, i){
  main = paste0(out$Symbol[i], ": ", out$gencodeID[i])
  p = paste0("p=",signif(out$P_PrimaryDxControl[i],3))
  titles <- c(main, p)
  names(titles) <- c("main", "p")
  return(titles)
}

residExp_beeswarm <- function(expDx_table, titles){
  ggplot(expDx_table, aes(Dx, exp, color = Dx)) + 
    geom_beeswarm() +
    scale_color_manual(values = mdd_Dx_colors) +
    labs(y = "Residualized Expression",
         title = titles[["main"]],
         subtitle = titles[["p"]])
}


topGenes_residExp_beeswarm <- function(out, exprs, pd, pdf_name , n = 100){
  sigOrderMat = as.data.frame(apply(out[,c("P_PrimaryDxControl","P_PrimaryDxBipolar")], 2, function(x) order(x)[1:n]))
  ooL = sigOrderMat$"P_PrimaryDxControl"
  
  pdf(pdf_name,h=6,w=6)
  for(i in ooL) {
    expDx_table <- make_expDx_table(exprs, i , pd)
    exp_title <- exp_title(out, i)
    print(residExp_beeswarm(expDx_table, exp_title))
  }
  dev.off()
}


