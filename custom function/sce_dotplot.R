sce_dotplot <- function(object, # scobj
                        assay = NULL, # Name of assay to use, defaults to the active assay
                        idents = NULL, # The column in the metadata of scobj
                        markers = NULL, # A list of genes
                        levels = NULL, # A list in a specific order
                        text_x_size = 10, # The size of x-axis text (Default: 10)
                        text_y_size = 10, # The size of y-axis text (Default: 10)
                        text_x_color = 'black', # The color of x-axis text (Default: 'black')
                        text_y_color = 'black', # The color of y-axis text (Default: 'black')
                        dot_legend = NULL, # The legend name of dot
                        dot_color_low = '#0da9ce', # The lower color of dot (Default: '#0da9ce')
                        dot_color_mid = '#FFF6CA', # The middle color of dot (Default: '#FFF6CA')
                        dot_color_high = '#e74a32', # The higher color of dot (Default: '#e74a32')
                        title = NULL, # The title of dotplot
                        title_size = NULL # The size of title
                        ){
  
  suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)})
  
  Idents(object) <- idents
  df <- DotPlot(object, features = markers, assay = assay)$data

  colnames(df) <- c('avg.exp', 'Percentage', 'Marker', 'Celltype', 'Mean.expression')
  
  df$Celltype <- factor(df$Celltype, 
                        levels = levels)  
  
  p <- ggplot(df, aes(Marker, Celltype)) + 
       theme_bw() +
       geom_point(aes(size = Percentage, 
                      color = Mean.expression)) +
       scale_color_gradient2(low = dot_color_low, 
                             mid = dot_color_mid, 
                             high = dot_color_high,
                             name = dot_legend)+
       theme(axis.text.x = element_text(vjust = 1,
                                        hjust = 1,
                                        angle = 45,
                                        size = text_x_size,
                                        color = text_x_color),
             axis.text.y = element_text(size = text_y_size,
                                        color = text_y_color),
             plot.title = element_text(size = title_size, 
                                       hjust = 0.5),
             legend.position = 'right') + 
       coord_flip()+
       xlab('') + 
       ylab('') +
       ggtitle(title)
  
}
