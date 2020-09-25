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

i <- 18843
test <- make_expDx_table(gSaccExprs, i, pdSacc)

exp_title <- function(out, i){
  main = paste0(out$Symbol[i], "\n", out$gencodeID[i])
  p = paste0("p=",signif(out$P_PrimaryDxControl[i],3))
  titles <- c(main, p)
  names(titles) <- c("main", "p")
  return(titles)
}
test_title <- exp_title(outGene_sACC, i)

residExp_beeswarm <- function(expDx_table, titles){
  ggplot(expDx_table, aes(Dx, exp, color = Dx)) + 
    geom_beeswarm() +
    # scale_color_manual(mdd_colors) +
    labs(y = "Residualized Expression",
         title = titles[["main"]],
         subtitle = titles[["p"]])
}

test_bee <- residExp_beeswarm(test, test_title)
ggsave("Beeswarm_test.png", test_bee)
