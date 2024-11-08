sce_volcanoplot <- function(object, # A dataframe including gene, cluster, p_val, avg_log2FC and p_val_adj
                            levels = NULL, # A list that is arranged according to the order of cluster
                            back_col = "grey85", # The color of background (Default: "grey85")
                            back_width = 0.8, # The width of background (Default: 0.8)
                            back_alpha = 0.5, # The transparency of background (Default: 0.5)
                            jitter_size = 3, # The size of dot (Default: 3)
                            jitter_width = 0.4, # The width of dot (Default: 0.4)
                            jitter_alpha = 0.4, # The transparency of dot (Default: 0.4)
                            text_size = 4.5, # The size of cluster text (Default: 4.5)
                            text_col = "black", # The color of cluster text (Default: "black")
                            title = NULL, # The title of volcano plot
                            group_col =  c('#e74a32','#0da9ce'), # The color codes for up- and down-regulated genes
                            cluster_col = NULL, # The color codes for cluster
                            ptype = c('pvalue', 'adjpvalue'), 
                            pvalue_cutoff = 0.05, # The cutoff of pvalue (default 0.05)
                            adjpvalue_cutoff = 0.05, # The cutoff of adjpvalue (default 0.05)
                            log2FC_cutoff = 0.25 # The cutoff of log2FoldChange (default 0.25)
                            ){
  
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(RColorBrewer)
    library(grid)
    library(scales)})
  
  if (ptype == "pvalue"){
      up <- object %>%
            filter(p_val < pvalue_cutoff & avg_log2FC > log2FC_cutoff) 
  
      up$group <- 'Significant upregulation'
    
      down <- object %>% 
              filter(p_val < pvalue_cutoff & avg_log2FC < -log2FC_cutoff) 
    
      down$group <- 'Significant downregulation'
  
      df <- rbind(up, down)
  }
    
  if (ptype == "adjpvalue"){
      up <- object %>%
            filter(p_val_adj < adjpvalue_cutoff & avg_log2FC > log2FC_cutoff) 
    
      up$group <- 'Significant upregulation'
    
      down <- object %>% 
              filter(p_val_adj < adjpvalue_cutoff & avg_log2FC < -log2FC_cutoff) 
    
      down$group <- 'Significant downregulation'
    
      df <- rbind(up, down)
  }
  
  df$group <- factor(df$group, levels = c('Significant upregulation', 'Significant downregulation'))

  df_bg <- df%>%
           group_by(cluster)%>%
           summarize(max_log2FC = max(avg_log2FC), 
                     min_log2FC = min(avg_log2FC)
                    )
  
  levels <- levels
  df_bg$cluster <- factor(df_bg$cluster, 
                          levels = levels)
  
  p <- ggplot()+
       geom_col(data = df_bg,
             mapping = aes(cluster, max_log2FC),
             fill = back_col,
             width = back_width,
             alpha = back_alpha) +
    
       geom_col(data = df_bg,
                mapping = aes(cluster, min_log2FC),
                fill = back_col,
                width = back_width,
                alpha = back_alpha) +
    
       geom_jitter(data = df,
                   mapping = aes(x = cluster, y = avg_log2FC, color = group),
                   size = jitter_size,
                   width = jitter_width,
                   alpha = jitter_alpha) +
    
       geom_col(data = df_bg,
                mapping = aes(x = cluster, y = 0.3, fill = cluster)) +
    
       geom_col(data = df_bg,
                mapping = aes(x = cluster, y = -0.3, fill = cluster)) +
    
       geom_text(data = df_bg,
                 mapping = aes(x = cluster, y = 0, label = cluster),
                 size = text_size,
                 color = text_col) +
       scale_color_manual(values = group_col)+
       scale_fill_manual(values = cluster_col) +
    
       theme_classic() +
       theme(axis.line.x = element_blank(),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             axis.line.y = element_line(linewidth = 0.8),
             axis.text.y = element_text(size = 12, color = "black"),
             axis.title.y = element_text(size = 20, face = 'bold', color = "black"),
             axis.ticks.y = element_line(linewidth = 0.8),
             plot.title = element_text(vjust = 0.5, hjust = 0.5, size = 20, face = 'bold'),
             legend.position = c(0.85, 1)) +
    labs(x = "", y = "Average Log2FC", fill = NULL, color = NULL)+
    ggtitle(title) +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1),
                                label.theme = element_text(size = 12)),
           fill = 'none')
  
  return(p)
}
