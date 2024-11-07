RunQC <- function(object,
                  org = 'hsa', 
                  LowerFeatureCutoff = 200,
                  UpperFeatureCutoff = "MAD",
                  UpperMitoCutoff = 10,
                  Hb = T, # percent.Hb
                  HbCutoff = 0,
                  decontX = T,
                  decontXCutoff = 0.2,
                  doubletdetection = T,
                  dir
){
  suppressPackageStartupMessages({
    library(Seurat)
    library(scds)
    library(dplyr)
    library(decontX)})
  if(UpperFeatureCutoff != "MAD" & !is.numeric(UpperFeatureCutoff)) {
    stop("Please use MAD and numeric cutoff for UpperFeatureCount!")
  }
  
  if(doubletdetection == T){
    sce= as.SingleCellExperiment(object)
    sce = cxds(sce, estNdbl = T)
    sce = bcds(sce, estNdbl = T)
    sce = cxds_bcds_hybrid(sce, estNdbl = T)
    object = as.Seurat(sce)
  }
  
  if(org == 'hsa'){
    object[["percent.mito"]] <- PercentageFeatureSet(object, pattern = "^MT-")
    object[["percent.rb"]] <- PercentageFeatureSet(object, pattern = "^RP[SL]")
    Hb_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    Hb_m <- match(Hb_genes, rownames(object@assays$RNA)) 
    Hb_genes <- rownames(object@assays$RNA)[Hb_m]
    Hb_genes <- Hb_genes[!is.na(Hb_genes)]  
    object[["percent.hb"]] <- PercentageFeatureSet(object, features = Hb_genes) 
    
  }else{
    object[["percent.mito"]] <- PercentageFeatureSet(object, pattern = "^mt-")
    object[["percent.rb"]] <- PercentageFeatureSet(object, pattern = "^Rp[sl]")
    Hb_genes <-  c('Hba-a1','Hba-a2','Hbb-bs','Hbb-bt','Hbb-b1')
    Hb_m <- match(Hb_genes, rownames(object@assays$RNA)) 
    Hb_genes <- rownames(object@assays$RNA)[Hb_m]
    Hb_genes <- Hb_genes[!is.na(Hb_genes)]  
    object[["percent.hb"]] <- PercentageFeatureSet(object, features = Hb_genes)
  }
  
  write.csv(object@meta.data, file=paste(dir,"/percentmito_hb_rb.csv",sep = ""))
  png(paste(dir, "/QC_Vlnplot.png",sep = ""), width = 10, height = 6, units = "in", res = 600)
  print({VlnPlot(object = object, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.hb", 'percent.rb'),
                 ncol = 5)})
  dev.off()
  
  object@misc[["filterstats"]] <- list()
  object@misc[["filterstats"]][['TotalCellsbeforefilteration']] <- dim(object)[2]
  

    #Using a median + 3 MAD cutoff for high genes.
    if(UpperFeatureCutoff == "MAD"){
      UpperFeatureCutoff <- median(object$nFeature_RNA) + 3*mad(object$nFeature_RNA)
      }
    #Detect percent.mito, percent.hb, nFeature_RNA, Doublet, ambient RNA
    object@misc[["filterstats"]][['TotalSamples']] <- dim(object[[]][1])[1]
    
    cells.use <- colnames(object)[which(object[[]]['percent.mito'] < UpperMitoCutoff)]
    object@misc[["filterstats"]][['Mitofilter']] <- dim(object[[]][1])[1] - length(cells.use)
    object <- subset(object, cells = cells.use)
    
    if(Hb == T){
    cells.use <- colnames(object)[which(object[[]]['percent.hb'] <= HbCutoff)]
    object@misc[["filterstats"]][['Hbfilter']] <- dim(object[[]][1])[1] - length(cells.use)
    object <- subset(object, cells = cells.use)
    }
    
    cells.use <- colnames(object)[which(object[[]]['nFeature_RNA'] > LowerFeatureCutoff)]
    object@misc[["filterstats"]][['LowFeatureFilter']] <- dim(object[[]][1])[1] - length(cells.use)
    object <- subset(object, cells = cells.use)
    
    cells.use <- colnames(object)[which(object[[]]['nFeature_RNA'] < UpperFeatureCutoff)]
    object@misc[["filterstats"]][['HighFeatureFilter']] <- dim(object[[]][1])[1] - length(cells.use)
    object <- subset(object, cells = cells.use)
    
    if(doubletdetection == T){
    cells.use <- rownames(object@meta.data %>% dplyr::filter(cxds_call == TRUE &
                                                             bcds_call == TRUE &
                                                             hybrid_call == TRUE))
    object@misc[["filterstats"]][['DoubletFilter']] <- length(cells.use)
    object <- subset(object, cells = setdiff(Cells(object), cells.use))
    }
    
    if(decontX == T){
    result <- decontX(object@assays$RNA@counts)
    object$contamination <- result$contamination
    cells.use <- colnames(object)[which(object[[]]['contamination'] < decontXCutoff)]
    object@misc[["filterstats"]][['DecontXFilter']] <- dim(object[[]][1])[1] - length(cells.use)
    object <- subset(object, cells = cells.use)
    }
   
    return(object)
}
