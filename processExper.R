processExper <- function(object, 
                         org = c("hsa", "mmu"),
                         cyclescoring = T,
                         sct.method = F,
                         reduction = "cca"){
  suppressPackageStartupMessages({
    library(Seurat)
    library(glmGamPoi)
    library(dplyr)
    library(tidyverse)})
  
  load('./mouseortholog.RData')
  m2h = mouse_human
  
  if(org != "hsa" & org != "mmu") {
    stop("Please input hsa or mmu for org!")
  }
  
  # Single scobj
  if(length(table(object$Sample)) < 2) {
    object <- object %>%
              NormalizeData(normalization.method = "LogNormalize",
                            scale.factor = 10000,
                            verbose = T) %>%
              FindVariableFeatures(selection.method = "vst", 
                                   nfeatures = 2000, 
                                   verbose = T) 
    if(cyclescoring == T){
      if(org == 'hsa'){
        # Assign scores in the CellCycleScoring function. Stores S and G2/M scores in object@meta.data, along with the predicted classification of each cell in either G2M or S phase
        object <- CellCycleScoring(object = object, 
                                   s.features = cc.genes.updated.2019$s.genes, 
                                   g2m.features = cc.genes.updated.2019$g2m.genes)
      }else{
        cc.genes.updated.2019$s.genes <- m2h %>% filter(human_name %in% cc.genes.updated.2019$s.genes) %>% pull(mouse_name)
        cc.genes.updated.2019$g2m.genes <- m2h %>% filter(human_name %in% cc.genes.updated.2019$g2m.genes) %>% pull(mouse_name)
        object <- CellCycleScoring(object = object, 
                                   s.features  = cc.genes.updated.2019$s.genes, 
                                   g2m.features = cc.genes.updated.2019$g2m.genes)
      }
    }
    # NormalizeData-normalize data
    # FindVariableFeatures-detection of variable genes; calculates the average expression and dispersion for each gene, places these genes into bins, and then calculates a z-score for dispersion within each bin
    # Scaling the data and removing unwanted sources of variation
    if(sct.method == T){
      if(cyclescoring == T){
        object <- SCTransform(object,
                              method = "glmGamPoi",
                              vars.to.regress = c("percent.mito", "S.Score", "G2M.Score"),
                              verbose = TRUE)
      }else{ 
        object <- SCTransform(object,
                              method = "glmGamPoi",
                              vars.to.regress = c("percent.mito"),
                              verbose = TRUE) 
      }
    }  
    
    if(sct.method == F){
      if(cyclescoring == T){
        object <- SCTransform(object,
                              vars.to.regress = c("percent.mito", "S.Score", "G2M.Score"),
                              verbose = TRUE)
      }else{ 
        object <- SCTransform(object,
                              vars.to.regress = c("percent.mito"),
                              verbose = TRUE) 
      }
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
      if(cyclescoring == T){
        if(org == 'hsa'){
          # Assign scores in the CellCycleScoring function.Stores S and G2/M scores in object@meta.data, along with the predicted classification of each cell in either G2M, S or G1 phase
          object[[i]] <- CellCycleScoring(object = object[[i]], 
                                          s.features = cc.genes.updated.2019$s.genes, 
                                          g2m.features = cc.genes.updated.2019$g2m.genes)
        }else{
          cc.genes.updated.2019$s.genes <- m2h %>% filter(human_name %in% cc.genes.updated.2019$s.genes) %>% pull(mouse_name)
          cc.genes.updated.2019$g2m.genes <- m2h %>% filter(human_name %in% cc.genes.updated.2019$g2m.genes) %>% pull(mouse_name)
          object[[i]] <- CellCycleScoring(object = object[[i]], 
                                          s.features  = cc.genes.updated.2019$s.genes, 
                                          g2m.features = cc.genes.updated.2019$g2m.genes)
        }
      }
    }
    
    for (i in 1:length(object)) {
      if(sct.method == T){
        if(cyclescoring == T){
          object[[i]] <- SCTransform(object[[i]], 
                                     method = "glmGamPoi",
                                     vars.to.regress = c("percent.mito", "S.Score", "G2M.Score"),
                                     verbose = TRUE)
        }else{
          object[[i]] <- SCTransform(object[[i]], 
                                     method = "glmGamPoi",
                                     vars.to.regress = c("percent.mito", "S.Score", "G2M.Score"),
                                     verbose = TRUE)
        }
      }
      
      if(sct.method == F){
        if(cyclescoring == T){
          object[[i]] <- SCTransform(object[[i]], 
                                     vars.to.regress = c("percent.mito", "S.Score", "G2M.Score"),
                                     verbose = TRUE)
        }else{
          object[[i]] <- SCTransform(object[[i]], 
                                     vars.to.regress = c("percent.mito"),
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

