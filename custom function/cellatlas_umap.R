cellatlas_umap <- function(object, # scobj
                           idents = NULL, # The column in the metadata of scobj
                           levels = NULL, # A list in a specific order
                           hull_alpha = NULL, # The transparency of hull
                           hull_size = NULL, # The size of hull curve
                           hull_lty = NULL, # The type of hull curve
                           hull_delta = NULL, # The distance to extend the hull curve
                           dot_color = NULL,  # A list of color codes
                           dot_size = NULL, # The size of dot
                           dot_alpha = NULL, # The transparency of dot
                           label_size = NULL, # The size of label
                           label_color = F # The color of label (Default: black)
                          ){
  
  suppressPackageStartupMessages({
    library(Seurat)
    library(ggrepel)
    library(ggunchull)
    library(ggplot2)})
  
  Idents(object) <- idents
  object@meta.data$Custom <- Idents(object)
  
  df <- object@reductions$umap@cell.embeddings %>%
        as.data.frame() %>% 
        cbind(cell_type = object@meta.data$Custom)
  
  df$cell_type <- factor(x = df$cell_type, 
                         levels = levels)
  
  colors <- dot_color
  
  p <- ggplot(df, aes(x = UMAP_1, 
                      y = UMAP_2, 
                      color = cell_type)) +
       stat_unchull(alpha = hull_alpha, 
                    size = hull_size, 
                    lty = hull_lty, 
                    delta = hull_delta) +
       geom_point(size = dot_size, 
                  alpha = dot_alpha)+
       scale_color_manual(values = colors)+
       theme(panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(), 
             panel.border = element_blank(), 
             axis.title = element_blank(),  
             axis.text = element_blank(), 
             axis.ticks = element_blank(),
             panel.background = element_rect(fill = 'white'), 
             plot.background = element_rect(fill="white"),
             legend.position = 'none' ) +
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
  
  med <- df %>%
         group_by(cell_type) %>%
         summarise(UMAP_1 = median(UMAP_1),
                   UMAP_2 = median(UMAP_2))
  
  if(label_color == T){
     p + geom_label_repel(aes(label = cell_type), 
                       data = med,
                       size = label_size) 
  }else{
     p + geom_label_repel(aes(label = cell_type), 
                         data = med,
                         size = label_size,
                         color = 'black')  
  }
}
