my_GSEAbarplot <- function(object, # A dataframe containing Pathway1, Pathway2 and NES columns
                           levels = NULL, # A list in a specific order
                           maxPathway = NULL,
                           title = NULL,
                           group_color = NULL,
                           pathway_color = NULL
                           ){
  
  suppressPackageStartupMessages({
    library(tidyverse)
    library(dplyr)
    library(forcats)
    library(ggplot2)})
  
  levels <- levels
  df <- object
  
  # Add Pathway column to Pathway2 column and set NES to 0
  df %>% distinct(Pathway1, .keep_all = T) %>%  mutate(Pathway2 = Pathway1, NES = 0) -> a
  
  # Merge, sort by Pathway and NES, and select the desired columns
  df <- rbind(df, a) %>% arrange(Pathway1, desc(NES)) %>% select(Pathway1, Pathway2, NES)
  
  # Add a Label column and mark the data as null for empty NES values
  df$Label <- ifelse(df$NES == 0, "", df$NES)
  
  # Sort by preselected order
  df$Pathway1 <- factor(df$Pathway1, levels = levels)
  df <- df[order(df$Pathway1), ]
  
  # Set levels
  df$Pathway2 <- df$Pathway2 %>% factor() %>% fct_inorder() %>% fct_rev()

  
  row_names <- df$Pathway2
  
  maxPathway <- maxPathway


  if(length(levels) == 1){
    row_indices1 <- which(row_names %in% levels[1])
    row_indices2 <- which(row_names %in% maxPathway[1])
    
    row1 <- df %>% slice(row_indices1)
    row2 <- df %>% slice(row_indices2)
    
    df1 <- df %>%
           slice(-1:-length(row_names)) %>%
           bind_rows(row1)  %>%  
           bind_rows(row2)  %>% 
           bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices1, -row_indices2))
    df <- df1
    }
  
  if(length(levels) == 2){
    row_indices1 <- which(row_names %in% levels[1])
    row_indices2 <- which(row_names %in% maxPathway[1])
    
    row1 <- df %>% slice(row_indices1)
    row2 <- df %>% slice(row_indices2)
    
    df1 <- df %>%
           slice(-1:-length(row_names)) %>%
           bind_rows(row1)  %>%  
           bind_rows(row2)  %>% 
           bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices1, -row_indices2))
    
    row_indices3 <- which(row_names %in% levels[2])
    row_indices4 <- which(row_names %in% maxPathway[2])
 
    row3 <- df %>% slice(row_indices3)
    row4 <- df %>% slice(row_indices4)
    
    df2 <- df1 %>%
      slice(-row_indices4:-length(row_names)) %>%
      bind_rows(row3)  %>%  # 先添加第 5 行，再添加第 3 行
      bind_rows(row4) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices3, -row_indices4) %>% slice(-1:-(row_indices4-1)))
    
    df <- df2
    } 
  
  if(length(levels) == 3){
    row_indices1 <- which(row_names %in% levels[1])
    row_indices2 <- which(row_names %in% maxPathway[1])
    
    row1 <- df %>% slice(row_indices1)
    row2 <- df %>% slice(row_indices2)
    
    df1 <- df %>%
      slice(-1:-length(row_names)) %>%
      bind_rows(row1)  %>%  
      bind_rows(row2)  %>% 
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices1, -row_indices2))
    
    row_indices3 <- which(row_names %in% levels[2])
    row_indices4 <- which(row_names %in% maxPathway[2])
    
    row3 <- df %>% slice(row_indices3)
    row4 <- df %>% slice(row_indices4)
    
    df2 <- df1 %>%
      slice(-row_indices4:-length(row_names)) %>%
      bind_rows(row3)  %>%  # 先添加第 5 行，再添加第 3 行
      bind_rows(row4) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices3, -row_indices4) %>% slice(-1:-(row_indices4-1)))
    
    row_indices5 <- which(row_names %in% levels[3])
    row_indices6 <- which(row_names %in% maxPathway[3])
    
    row5 <- df %>% slice(row_indices5)
    row6 <- df %>% slice(row_indices6)
    
    df3 <- df2 %>%
      slice(-row_indices6:-length(row_names)) %>%
      bind_rows(row5)  %>%  # 先添加第 5 行，再添加第 3 行
      bind_rows(row6) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices5, -row_indices6) %>% slice(-1:-(row_indices6-1)))
    
    df <- df3
  } 
  
  
  if(length(levels) == 4){
    row_indices1 <- which(row_names %in% levels[1])
    row_indices2 <- which(row_names %in% maxPathway[1])
    
    row1 <- df %>% slice(row_indices1)
    row2 <- df %>% slice(row_indices2)
    
    df1 <- df %>%
      slice(-1:-length(row_names)) %>%
      bind_rows(row1)  %>%  
      bind_rows(row2)  %>% 
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices1, -row_indices2))
    
    row_indices3 <- which(row_names %in% levels[2])
    row_indices4 <- which(row_names %in% maxPathway[2])
    
    row3 <- df %>% slice(row_indices3)
    row4 <- df %>% slice(row_indices4)
    
    df2 <- df1 %>%
      slice(-row_indices4:-length(row_names)) %>%
      bind_rows(row3)  %>%  # 先添加第 5 行，再添加第 3 行
      bind_rows(row4) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices3, -row_indices4) %>% slice(-1:-(row_indices4-1)))
    
    row_indices5 <- which(row_names %in% levels[3])
    row_indices6 <- which(row_names %in% maxPathway[3])
    
    row5 <- df %>% slice(row_indices5)
    row6 <- df %>% slice(row_indices6)
    
    df3 <- df2 %>%
      slice(-row_indices5:-length(row_names)) %>%
      bind_rows(row5)  %>%  # 先添加第 5 行，再添加第 3 行
      bind_rows(row6) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices5, -row_indices6) %>% slice(-1:-(row_indices6-1)))
    
    row_indices7 <- which(row_names %in% levels[4])
    row_indices8 <- which(row_names %in% maxPathway[4])
    
    row7 <- df %>% slice(row_indices7)
    row8 <- df %>% slice(row_indices8)
    
    df4 <- df3 %>%
      slice(-row_indices7:-length(row_names)) %>%
      bind_rows(row7)  %>%  # 先添加第 5 行，再添加第 3 行
      bind_rows(row8) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices7, -row_indices8) %>% slice(-1:-(row_indices8-1)))
    
    df <- df4
  } 
  
  
  
  df$Pathway2 <- fct_inorder(df$Pathway2) %>% fct_rev()
  
  barplot <- ggplot(df, aes(x = NES, Pathway2, fill = Pathway1))+
             geom_col()+ 
             geom_text(aes(label = Label), 
                       color = 'black',
                       size = 2.5,
                       hjust = "left",
                       nudge_x = 0.1) + 
             scale_x_continuous(limits = c(0, 2.3),
                                expand = expansion(mult = c(0, .1))) + 
             labs(x = "NES", y = "", title = title) +
             scale_fill_manual(values = pathway_color)
  
  g <- ggplot_build(barplot)
  mycol <- g$data[[1]]["fill"]
  col <- rev(mycol$fill)

  num <- rev(df$NES)
  index <- which(num == 0)
  col[index] <- group_color
  

  theme <-  theme_bw()+
            theme(plot.title = element_text(size = rel(1),
                                            hjust = 0.5,
                                            face = 'plain'),
                  axis.title = element_text(size = rel(1)),
                  axis.text.x = element_text(color = 'black',size=rel(1)),
                  axis.text.y = element_text(size = rel(1),
                                             color = col,
                                             face = 'bold'),
                  legend.position = "none",
                  plot.margin = unit(x = c(top.mar = 0.2,
                                           right.mar = 0.2,
                                           left.mar=0.2,
                                           bottom.mar= 0.2),
                                     units = 'inches'))
  
  barplot1 <- barplot + theme
}
