FindMarker_genes <- function(dataset, 
                             clusters, 
                             comparison,
                             ...) {

  suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(dplyr)})

  Idents(dataset) <- paste(as.character(Idents(dataset)), dataset[[comparison[1]]][[1]], sep = "_")
  
  de <- list()
  
  for (i in seq_along(1:length(clusters))) {
    
    d <- FindMarkers(dataset, 
                     ident.1 = paste(clusters[i], comparison[3], sep = "_"),
                     ident.2 = paste(clusters[i], comparison[2], sep = "_"),
                     assay = "SCT",
                     ...)
    de[[i]] <- d %>%
      rownames_to_column(var = "gene") %>%
      add_column(cluster = clusters[i], .after = 1) 
    
    
  }
  return(do.call("rbind", de))
}
