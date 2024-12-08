Split_Vln_stacked <- function(object, # scobj
                              assay, # Name of assay to use, defaults to the active assay
                              feature, # A list of genes
                              split.by = NULL, # Split clusters into different groups
                              split.plot = F, # Split VlnPlot
                              pt.size, # The size of dot
                              face = 'bold', # # The font of x-axis text and y-axis title (Default: bold)
                              text_x_size = 10, # The size of x-axis text (Default: 10)
                              text_y_size = 10, # The size of y-axis text (Default: 10)
                              title_y_size = 10, # The size of y-axis title (Default: 10)
                              cols # The color of clusters
                              )
  {
  suppressPackageStartupMessages({
    library(Seurat)
    library(dittoSeq)
    library(deeptime)
    library(ggplot2)
    library(ggpubr)})
  
  p1list <- list()
  for (i in 1:length(feature)){
      p1 <- VlnPlot(object,  
                    assay = assay,
                    features = feature[i], 
                    split.by = split.by,
                    split.plot = split.plot, 
                    pt.size = pt.size)+
        theme_bw()+
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = text_y_size),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size = title_y_size,
                                          colour = 'black',
                                          face = face,
                                          angle = 0,
                                          vjust = 0.5),
              axis.ticks.x = element_blank(),
              title = element_blank(),
              legend.position = 'none',
              panel.border = element_rect(fill = NA),
              plot.margin = margin(-0.05, 0, -0.05, 0, "cm"),
              panel.grid = element_blank())+
        ylab(feature[i])+
        scale_fill_manual(values = cols)
      p1list[[i]] <- p1
  }
  
  if(split.plot == T){
    p2 <- VlnPlot(object,  
                  assay = assay, 
                  features = tail(feature, 1), 
                  split.by = split.by,
                  split.plot = split.plot, 
                  pt.size = pt.size)+
      theme_bw()+
      theme(axis.text.y = element_text(size = text_y_size),
            axis.text.x = element_text(size = text_x_size, 
                                       colour = 'black',
                                       face = face,
                                       angle = 45, 
                                       vjust = 1, 
                                       hjust = 1),
            axis.title.y = element_text(size = title_y_size,
                                        colour = 'black',
                                        face = face, 
                                        angle = 0, 
                                        vjust = 0.5),
            axis.ticks = element_blank(),
            title = element_blank(),
            legend.position = 'bottom',
            panel.border = element_rect(fill = NA),
            plot.margin = margin(-0.05, 0, 0, 0, "cm"),
            panel.grid = element_blank(),
            legend.box.background = element_blank(),
            legend.text = element_text(color = "black",
                                       size=10),
            legend.spacing.x = unit(0.2,'cm'),
            legend.key.width = unit(0.4,'cm'),
            legend.key.height = unit(0.4,'cm'),
            legend.background = element_blank())+
      ylab(tail(feature, 1))+
      scale_fill_manual(values = cols)
    p22 <- p2 + ylim(0, max(p2$data[tail(feature,1)] + 1.5))
  }else{
    p2 <- VlnPlot(object, 
                  assay = assay,
                  features = tail(feature, 1),
                  split.by = split.by,
                  split.plot = split.plot,
                  pt.size = pt.size)+
      theme_bw()+
      theme(axis.text.y = element_text(size = text_y_size),
            axis.text.x = element_text(size = text_x_size, 
                                       colour = 'black',
                                       face = face,
                                       angle = 45, 
                                       vjust = 1, 
                                       hjust = 1),
            axis.title.y = element_text(size = title_y_size, 
                                        colour = 'black',
                                        face = face,
                                        angle = 0, 
                                        vjust = 0.5),
            axis.ticks.x = element_blank(),
            title = element_blank(),
            legend.position = 'none',
            panel.border = element_rect(fill = NA),
            plot.margin = margin(-0.05, 0, 0, 0, "cm"),
            panel.grid = element_blank())+
      ylab(tail(feature, 1))+
      scale_fill_manual(values = cols)
    p22 <- p2 + ylim(0, max(p2$data[tail(feature,1)] + 1.5))
  }
    p1list[[length(feature)]] <- p2
    
  p <- deeptime::ggarrange2(plots = p1list, nrow = length(feature))
  return(p)
}
