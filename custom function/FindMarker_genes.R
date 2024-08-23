FindMarker_genes <- function(object, # scobj
                             assay = NULL, # Name of assay to use, defaults to the active assay
                             clusters, # A list of preselected clusters
                             comparison, # Format: c('Condition', 'Control', 'Case')
                             ...) {

  suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(dplyr)})

  Idents(object) <- paste(as.character(Idents(object)), object[[comparison[1]]][[1]], sep = "_")
  
  de <- list()
  
  for (i in seq_along(1:length(clusters))) {
    
    d <- FindMarkers(object, 
                     ident.1 = paste(clusters[i], comparison[3], sep = "_"),
                     ident.2 = paste(clusters[i], comparison[2], sep = "_"),
                     assay = assay,
                     ...)
    de[[i]] <- d %>%
               rownames_to_column(var = "gene") %>%
               add_column(cluster = clusters[i], .after = 1) 
    
    
  }
  return(do.call("rbind", de))
}
