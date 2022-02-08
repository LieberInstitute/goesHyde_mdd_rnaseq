
get_signif <- function(outFeature, colname = "common_feature_id", cutoff = 0.05, return_unique = FALSE){
  signif <- outFeature[[colname]][outFeature$adj.P.Val < cutoff]
  if(return_unique) signif <- unique(signif)
  signif <- signif[!is.na(signif)]
  return(signif)
}

my_flatten <- function (x, use.names = TRUE, classes = "ANY") 
{
  #' Source taken from rlist::list.flatten
  len <- sum(rapply(x, function(x) 1L, classes = classes))
  y <- vector("list", len)
  i <- 0L
  items <- rapply(x, function(x) {
    i <<- i + 1L
    y[[i]] <<- x
    TRUE
  }, classes = classes)
  if (use.names && !is.null(nm <- names(items))) 
    names(y) <- nm
  return(y)
}


compare_tstat_plot <- function(DE_out_x, DE_out_y, name_x, name_y, title = "", save_png = FALSE, prefix = "", fdr_cutoff = 0.05){
  
  if(!identical(colnames(DE_out_x), colnames(DE_out_y)))  message("colnames don't match")
  if(!identical(rownames(DE_out_x), rownames(DE_out_y)))  message("rownames don't match")
  
  colnames(DE_out_x) <- paste0("x.", colnames(DE_out_x))
  colnames(DE_out_y) <- paste0("y.", colnames(DE_out_y))
  
  DE_out <- cbind(DE_out_x, DE_out_y)
  rho <- cor(DE_out$x.t, DE_out$y.t, method = "spearman")
  
  rho_anno = paste0("rho = ", format(round(rho, 2), nsmall = 2))
  # message(rho_anno)
  
 DE_out <- DE_out %>%
    mutate(Signif = case_when(x.adj.P.Val < fdr_cutoff & y.adj.P.Val < fdr_cutoff ~"sig_Both",
                              x.adj.P.Val < fdr_cutoff ~ paste0("sig_", name_x),
                              y.adj.P.Val < fdr_cutoff ~ paste0("sig_", name_y),
                              TRUE ~ "None"))
 
 sig_colors <- c(RColorBrewer::brewer.pal(3, "Set1"),"dark grey")
 names(sig_colors) <- c("sig_Both", paste0("sig_", name_x), paste0("sig_", name_y), "None")
 
 if(title == "") title = paste(name_x, "vs.", name_y)
 
 t_stat_plot <- ggplot(DE_out, aes(x = x.t, y = y.t, color = Signif)) +
   geom_point(alpha = 0.7, size = 0.5) +
   labs(x = paste("t-stats", name_x), 
        y = paste("t-stats", name_y),
        title = title, 
        subtitle = rho_anno, 
        parse = T) +
   scale_color_manual(values = sig_colors)+
   theme_bw()
 
 message(title)
 return(t_stat_plot)
  # 
  # 
  # if(save_png){
  #   fn = here("differential_expression","plots",paste0(prefix,".png"))
  #   ggsave(model_plot, filename = fn)              
  # }
  # print(model_plot)
}
