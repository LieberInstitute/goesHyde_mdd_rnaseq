
findMarkers_1vAll <- function(sce, assay = "logcounts", cellType_col = "cellType"){
  ct <- levels(sce[[cellType_col]])
  names(ct) <- ct
  ## Traditional t-test with design as in PB'd/limma approach ===
  mod <- with(colData(sce), model.matrix(~ donor))
  mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`
  
  markers.t.1vAll <- map(ct, function(x){
    sce$contrast <- ifelse(sce[[cellType_col]]==x, 1, 0)
    fm <- findMarkers(sce, groups=sce$contrast,
                assay.type=assay, design=mod, test="t",
                direction="up", pval.type="all", full.stats=T)
    fm <- fm[[2]]$stats.0
    fm.std <- findMarkers(sce, groups=sce$contrast,
                                     assay.type="logcounts", design=mod, test="t",
                                     std.lfc=TRUE,
                                     direction="up", pval.type="all", full.stats=T)
    fm.std <- fm.std[[2]]$stats.0
    colnames(fm.std)[[1]] <- "std.logFC"
    return(cbind(fm, fm.std[,1, drop=FALSE]))
  })
  
  markers.t.1vAll.table <- do.call("rbind",markers.t.1vAll) %>% 
    as.data.frame() %>%
    rownames_to_column("gene")%>% 
    as_tibble %>%
    add_column(cellType.target = rep(ct, each = nrow(sce)))
  
  return(markers.t.1vAll.table)
}
