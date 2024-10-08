sce_GSEA <- function(object, # A dataframe of DEGs in each cell cluster
                     clusters = NULL, # A list of preselected clusters
                     pathway # MSigDB pathways
                    ){

  suppressPackageStartupMessages({
    library(msigdbr)
    library(fgsea)
    library(tibble)
    library(tidyverse)
    library(dplyr)})

  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop(paste("Package \"fgsea\" needed for this function to work. Please install it."),
         call. = FALSE)
  }
  
  if (is.null(clusters)) clusters <- unique(object$cluster)
  
  gsea_res <- list()
  
  for (i in seq(length(clusters))) {
    
    data <- object %>%
            dplyr::filter(.data$cluster == clusters[i]) %>%
            arrange(desc(.data$avg_log2FC))
    
    l <- data$avg_log2FC
    names(l) <- data$gene
    
    res <- fgsea::fgsea(pathways = pathway, 
                        stats = l, 
                        minSize = 5, 
                        maxSize = 500, 
                        nperm = 10000)
    
    res <- res %>%
           add_column(cluster = clusters[i],
                      .before = 1)
    
    gsea_res[[i]] <- res
  }
  return(as_tibble(do.call("rbind", gsea_res)))
}
