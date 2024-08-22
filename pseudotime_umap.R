pseudotime_umap <- function(object, # scobj
                            pseudotime_obj = NULL, # Pseudotime obj
                            idents = NULL,  # The column in the metadata of scobj
                            sub_idents = NULL, # The pseudotime column in the metadata of scobj (eg. slingPseudotime_1)
                            dot_color = c("#440154FF", "#FDE725FF"), # A list of color codes
                            dot_size = NULL, # The size of dot
                            dot_alpha = NULL, # The transparency of dot
                            label_size = NULL # The size of label
                            ){
  
  suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)})
  
  pseudotime_obj@colData %>% colnames() %>% grep(pattern = "^sling", x = ., value = TRUE) -> slingshot_res_colnames
  for (cn in slingshot_res_colnames) {
    message(cn)
    object[[cn]] <- pseudotime_obj@colData[, cn]
  }
  
  Idents(object) <- sub_idents
  object@meta.data$Custom <- Idents(object)
  
  df <- object@reductions$umap@cell.embeddings %>%
        as.data.frame() %>% 
        cbind(cell_type = object@meta.data$Custom)

  df$Pseudotime <- as.character(df$cell_type) %>% as.numeric()

  colors <- dot_color
  
  p <- ggplot(df, aes(x = UMAP_1, 
                      y = UMAP_2, 
                      color = Pseudotime)) +
    geom_point(size = dot_size, 
               alpha = dot_alpha) +
    scale_color_gradient(low = colors[1], 
                         high = colors[2]) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          axis.title = element_blank(),  
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = 'white'), #背景色
          plot.background = element_rect(fill="white"),
          legend.position = 'right' ) +
    geom_segment(aes(x = min(df$UMAP_1) , 
                     y = min(df$UMAP_2) ,
                     xend = min(df$UMAP_1) + 3, 
                     yend = min(df$UMAP_2)),
                 colour = "black", 
                 size=0.5, 
                 arrow = arrow(length = unit(0.3,"cm")))+ 
    geom_segment(aes(x = min(df$UMAP_1), 
                     y = min(df$UMAP_2)  ,
                     xend = min(df$UMAP_1), 
                     yend = min(df$UMAP_2) + 3),
                 colour = "black", 
                 size = 0.5,
                 arrow = arrow(length = unit(0.3,"cm"))) +
    annotate("text", 
             x = min(df$UMAP_1) + 1.5, 
             y = min(df$UMAP_2) -1, 
             label = "UMAP_1",
             color="black", 
             size = 3, 
             fontface = "plain" ) + 
    annotate("text", 
             x = min(df$UMAP_1) - 1, 
             y = min(df$UMAP_2) + 1.5, 
             label = "UMAP_2",
             color = "black", 
             size = 3, 
             fontface = "plain", 
             angle = 90) 
  
  
  Idents(object) <- idents
  object@meta.data$Custom1 <- Idents(object)
  
  df1 <- object@reductions$umap@cell.embeddings %>%
         as.data.frame() %>% 
         cbind(cell_type = object@meta.data$Custom1)
  
  med <- df1 %>%
         group_by(cell_type) %>%
         summarise(UMAP_1 = median(UMAP_1),
                   UMAP_2 = median(UMAP_2))
  
  p + geom_label_repel(aes(label = cell_type), 
                       data = med,
                       size = label_size,
                       color = 'black')  
  
}
