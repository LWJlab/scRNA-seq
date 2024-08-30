Split_Vln_stacked <- function(object, # scobj
                              assay, # Name of assay to use, defaults to the active assay
                              feature, # A list of genes
                              split.by = NULL, # Split clusters into different groups
                              split.plot = NULL, # Split VlnPlot
                              pt.size, # The size of dot
                              size, # The size of axis text
                              cols, # The color of clusters
                              test = F, # Statistical analysis
                              test_method = NULL, # c('t.test','wilcox.test')
                              sig_label= NULL) # Add * or p value, c("p.signif","p.format"))
{
  
  suppressPackageStartupMessages({
    library(Seurat)
    library(dittoSeq)
    library(ggplot2)
    library(ggpubr)})
  
  p1list <- list()
  for (i in 1:length(feature)){
    if(test==T){
      p1 <- VlnPlot(object, 
                    assay = assay, 
                    features = feature[i], 
                    split.by = split.by,
                    split.plot = split.plot, 
                    pt.size = pt.size)+
            theme_bw()+
            theme(axis.text = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(size = size, 
                                              angle = 0, 
                                              vjust = 0.5),
                  axis.ticks = element_blank(),
                  title = element_blank(),
                  legend.position = 'none',
                  panel.border = element_rect(fill = NA),
                  plot.margin = margin(-0.05, 0, -0.05, 0, "cm"),
                  panel.grid = element_blank())+
            ylab(feature[i])+
            scale_fill_manual(values = cols)
      p11 <- p1 + ylim(0, max(p1$data[feature[i]] + 1.5))+
                  stat_compare_means(aes(group = split),
                                     label = sig_label,
                                     label.y = max(p1$data[feature[i]]),
                                     hide.ns = T)
      p1list[[i]] <- p11
    }else{
      p1 <- VlnPlot(object,  
                    assay = assay,
                    features = feature[i], 
                    split.by = split.by,
                    split.plot = split.plot, 
                    pt.size = pt.size)+
            theme_bw()+
            theme(axis.text = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(size = size,
                                              angle = 0,
                                              vjust = 0.5),
                  axis.ticks = element_blank(),
                  title = element_blank(),
                  legend.position = 'none',
                  panel.border = element_rect(fill = NA),
                  plot.margin = margin(-0.05, 0, -0.05, 0, "cm"),
                  panel.grid = element_blank())+
            ylab(feature[i])+
            scale_fill_manual(values = cols)
      p1list[[i]] <- p1
    }
  }
  if(split.plot==T){
    p2 <- VlnPlot(object,  assay = assay, features = tail(feature,1), split.by = split.by,
                  split.plot=split.plot, pt.size = pt.size)+
      theme_bw()+
      theme(axis.text.y = element_blank(),
            axis.text.x = element_text(colour = 'black', size = size, angle = 45, vjust = 1, hjust = 1),
            axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
            axis.ticks = element_blank(),
            title = element_blank(),
            legend.position = 'bottom',
            panel.border = element_rect(fill = NA),
            plot.margin = margin(-0.05,0,0,0, "cm"),
            panel.grid = element_blank(),
            legend.box.background = element_blank(),
            legend.text = element_text(color="black",size=10),
            legend.spacing.x=unit(0.2,'cm'),
            legend.key.width=unit(0.4,'cm'),
            legend.key.height=unit(0.4,'cm'),
            legend.background=element_blank())+
      ylab(tail(feature,1))+
      scale_fill_manual(values = cols)
    p22 <- p2+ylim(0, max(p2$data[tail(feature,1)]+1.5))+
      stat_compare_means(aes(group = split),
                         label =sig_label,
                         label.y = max(p2$data[tail(feature,1)]),
                         hide.ns=T)
  }else{
    p2 <- VlnPlot(object,  assay = assay, features = tail(feature,1), split.by = split.by,
                  split.plot=split.plot, pt.size = pt.size)+
      theme_bw()+
      theme(axis.text.y = element_blank(),
            axis.text.x = element_text(colour = 'black', size = size, angle = 45, vjust = 1, hjust = 1),
            axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
            axis.ticks = element_blank(),
            title = element_blank(),
            legend.position = 'none',
            panel.border = element_rect(fill = NA),
            plot.margin = margin(-0.05,0,0,0, "cm"),
            panel.grid = element_blank())+
      ylab(tail(feature,1))+
      scale_fill_manual(values = cols)
    p22 <- p2+ylim(0, max(p2$data[tail(feature,1)]+1.5))+
      stat_compare_means(aes(group = split),
                         label =sig_label,
                         label.y = max(p2$data[tail(feature,1)]),
                         hide.ns=T)
  }
  if(test==T){
    p1list[[length(feature)]] <- p22
  }else{
    p1list[[length(feature)]] <- p2
  }
  p<- deeptime::ggarrange2(plots = p1list, nrow = length(feature))
  return(p)
}
