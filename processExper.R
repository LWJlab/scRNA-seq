processExper <- function(object, 
                         sct.method = F,
                         reduction = "cca"){
  
  suppressPackageStartupMessages({
    library(Seurat)
    library(glmGamPoi)
    library(dplyr)
    library(tidyverse)})
  
  # Single scobj
  if(length(table(object$Sample)) < 2) {
    object <- object %>%
              NormalizeData(normalization.method = "LogNormalize",
                            scale.factor = 10000,
                            verbose = T) %>%
              FindVariableFeatures(selection.method = "vst", 
                                   nfeatures = 2000, 
                                   verbose = T) 
    
    # NormalizeData-normalize data
    # FindVariableFeatures-detection of variable genes; calculates the average expression and dispersion for each gene, places these genes into bins, and then calculates a z-score for dispersion within each bin
    # Scaling the data and removing unwanted sources of variation
    if(sct.method == T){
        object <- SCTransform(object,
                              method = "glmGamPoi",
                              vars.to.regress = c("percent.mito", "S.Score", "G2M.Score"),
                              verbose = TRUE)
      }else{ 
        object <- SCTransform(object,
                              vars.to.regress = c("percent.mito", "S.Score", "G2M.Score"),
                              verbose = TRUE)   
      }   
    
    # Multiple scobjs
  }else{ 
     object <- SplitObject(object, split.by = 'Sample')
    
    if(reduction != "cca" & reduction != "rpca" & reduction != "rlsi" ) {
      stop("Please input cca, rpca or rlsi for reduction!")
    }
    
    for (i in 1:length(object)) {
      object[[i]] <- NormalizeData(object[[i]],
                                   normalization.method = "LogNormalize",
                                   scale.factor = 10000,
                                   verbose = T) %>%
                     FindVariableFeatures(selection.method = "vst", 
                                          nfeatures = 2000, 
                                          verbose = T) 
    }
    
    for (i in 1:length(object)) {
      if(sct.method == T){
          object[[i]] <- SCTransform(object[[i]], 
                                     method = "glmGamPoi",
                                     vars.to.regress = c("percent.mito", "S.Score", "G2M.Score"),
                                     verbose = TRUE)
        }else{
          object[[i]] <- SCTransform(object[[i]], 
                                     vars.to.regress = c("percent.mito", "S.Score", "G2M.Score"),
                                     verbose = TRUE)
        } 
      }
    }
    
    
    features <- SelectIntegrationFeatures(object.list = object,
                                          nfeatures = 3000)
    
    object <- PrepSCTIntegration(object.list = object, 
                                 anchor.features = features)  
    
    anchors <- FindIntegrationAnchors(object.list = object,
                                      normalization.method = "SCT", 
                                      reduction = reduction,
                                      anchor.features = features)
    
    object <- IntegrateData(anchorset = anchors,
                            normalization.method = "SCT")
    
  } 
  object <- LogSeuratCommand(object = object)
}

