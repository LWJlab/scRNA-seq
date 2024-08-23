sce_GSEAbarplot <- function(object, # A dataframe containing Pathway1, Pathway2 and NES columns
                            levels = NULL, # A list that is arranged according to the order of Pathway1
                            maxPathway = NULL, # A list of the pathways with the highest NES values in each category within Pathway2.
                            category_color = 'grey10', # The color codes for category of pathways (Default: 'grey10')
                            pathway_color = NULL, # A list of pathways and their corresponding color codes
                            title = NULL, # The title of dotplot
                            title_size = 10, # The size of title (Default: 10)
                            num_size = 2.5, # The size of number (Default: 2.5)
                            text_x_size = 10, # The size of x-axis text (Default: 10)
                            text_y_size = 10, # The size of y-axis text (Default: 10)
                            xlim = c(0, 2.5) # A list of x-axis boundary
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
  
  maxPathway <- maxPathway
  
  if(length(levels) == 1){
    raw_row_names <- df$Pathway1
    raw_row <- levels[1]
    raw_row_indices <- which(raw_row_names %in% raw_row)
    
    df <- df %>% slice(raw_row_indices)
    
    row_names <- df$Pathway2
    
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
    raw_row_names <- df$Pathway1
    raw_row <- levels[1:2]
    raw_row_indices <- which(raw_row_names %in% raw_row)
    
    df <- df %>% slice(raw_row_indices)
    
    row_names <- df$Pathway2
    
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
      bind_rows(row3)  %>%  
      bind_rows(row4) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices3, -row_indices4) %>% slice(-1:-(row_indices4-1)))
    
    df <- df2
  } 
  
  if(length(levels) == 3){
    raw_row_names <- df$Pathway1
    raw_row <- levels[1:3]
    raw_row_indices <- which(raw_row_names %in% raw_row)
    
    df <- df %>% slice(raw_row_indices)
    
    row_names <- df$Pathway2
    
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
      bind_rows(row3)  %>%  
      bind_rows(row4) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices3, -row_indices4) %>% slice(-1:-(row_indices4-1)))
    
    row_indices5 <- which(row_names %in% levels[3])
    row_indices6 <- which(row_names %in% maxPathway[3])
    
    row5 <- df %>% slice(row_indices5)
    row6 <- df %>% slice(row_indices6)
    
    df3 <- df2 %>%
      slice(-row_indices6:-length(row_names)) %>%
      bind_rows(row5)  %>%  
      bind_rows(row6) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices5, -row_indices6) %>% slice(-1:-(row_indices6-1)))
    
    df <- df3
  } 
  
  
  if(length(levels) == 4){
    raw_row_names <- df$Pathway1
    raw_row <- levels[1:4]
    raw_row_indices <- which(raw_row_names %in% raw_row)
    
    df <- df %>% slice(raw_row_indices)
    
    row_names <- df$Pathway2
    
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
      bind_rows(row3)  %>%  
      bind_rows(row4) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices3, -row_indices4) %>% slice(-1:-(row_indices4-1)))
    
    row_indices5 <- which(row_names %in% levels[3])
    row_indices6 <- which(row_names %in% maxPathway[3])
    
    row5 <- df %>% slice(row_indices5)
    row6 <- df %>% slice(row_indices6)
    
    df3 <- df2 %>%
      slice(-row_indices6:-length(row_names)) %>%
      bind_rows(row5)  %>%  
      bind_rows(row6) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices5, -row_indices6) %>% slice(-1:-(row_indices6-1)))
    
    row_indices7 <- which(row_names %in% levels[4])
    row_indices8 <- which(row_names %in% maxPathway[4])
    
    row7 <- df %>% slice(row_indices7)
    row8 <- df %>% slice(row_indices8)
    
    df4 <- df3 %>%
      slice(-row_indices8:-length(row_names)) %>%
      bind_rows(row7)  %>%  
      bind_rows(row8) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices7, -row_indices8) %>% slice(-1:-(row_indices8-1)))
    
    df <- df4
  } 
  
  
  if(length(levels) == 5){
    raw_row_names <- df$Pathway1
    raw_row <- levels[1:5]
    raw_row_indices <- which(raw_row_names %in% raw_row)
    
    df <- df %>% slice(raw_row_indices)
    
    row_names <- df$Pathway2
    
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
      bind_rows(row3)  %>%  
      bind_rows(row4) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices3, -row_indices4) %>% slice(-1:-(row_indices4-1)))
    
    row_indices5 <- which(row_names %in% levels[3])
    row_indices6 <- which(row_names %in% maxPathway[3])
    
    row5 <- df %>% slice(row_indices5)
    row6 <- df %>% slice(row_indices6)
    
    df3 <- df2 %>%
      slice(-row_indices6:-length(row_names)) %>%
      bind_rows(row5)  %>%  
      bind_rows(row6) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices5, -row_indices6) %>% slice(-1:-(row_indices6-1)))
    
    row_indices7 <- which(row_names %in% levels[4])
    row_indices8 <- which(row_names %in% maxPathway[4])
    
    row7 <- df %>% slice(row_indices7)
    row8 <- df %>% slice(row_indices8)
    
    df4 <- df3 %>%
      slice(-row_indices8:-length(row_names)) %>%
      bind_rows(row7)  %>%  
      bind_rows(row8) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices7, -row_indices8) %>% slice(-1:-(row_indices8-1)))
    
    row_indices9 <- which(row_names %in% levels[5])
    row_indices10 <- which(row_names %in% maxPathway[5])
    
    row9 <- df %>% slice(row_indices9)
    row10 <- df %>% slice(row_indices10)
    
    df5 <- df4 %>%
      slice(-row_indices10:-length(row_names)) %>%
      bind_rows(row9)  %>%  
      bind_rows(row10) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices9, -row_indices10) %>% slice(-1:-(row_indices10-1)))
    
    df <- df5
  } 
  
  
  if(length(levels) == 6){
    raw_row_names <- df$Pathway1
    raw_row <- levels[1:6]
    raw_row_indices <- which(raw_row_names %in% raw_row)
    
    df <- df %>% slice(raw_row_indices)
    
    row_names <- df$Pathway2
    
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
      bind_rows(row3)  %>%  
      bind_rows(row4) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices3, -row_indices4) %>% slice(-1:-(row_indices4-1)))
    
    row_indices5 <- which(row_names %in% levels[3])
    row_indices6 <- which(row_names %in% maxPathway[3])
    
    row5 <- df %>% slice(row_indices5)
    row6 <- df %>% slice(row_indices6)
    
    df3 <- df2 %>%
      slice(-row_indices6:-length(row_names)) %>%
      bind_rows(row5)  %>% 
      bind_rows(row6) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices5, -row_indices6) %>% slice(-1:-(row_indices6-1)))
    
    row_indices7 <- which(row_names %in% levels[4])
    row_indices8 <- which(row_names %in% maxPathway[4])
    
    row7 <- df %>% slice(row_indices7)
    row8 <- df %>% slice(row_indices8)
    
    df4 <- df3 %>%
      slice(-row_indices8:-length(row_names)) %>%
      bind_rows(row7)  %>%  
      bind_rows(row8) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices7, -row_indices8) %>% slice(-1:-(row_indices8-1)))
    
    row_indices9 <- which(row_names %in% levels[5])
    row_indices10 <- which(row_names %in% maxPathway[5])
    
    row9 <- df %>% slice(row_indices9)
    row10 <- df %>% slice(row_indices10)
    
    df5 <- df4 %>%
      slice(-row_indices10:-length(row_names)) %>%
      bind_rows(row9)  %>%  
      bind_rows(row10) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices9, -row_indices10) %>% slice(-1:-(row_indices10-1)))
    
    row_indices11 <- which(row_names %in% levels[6])
    row_indices12 <- which(row_names %in% maxPathway[6])
    
    row11 <- df %>% slice(row_indices11)
    row12 <- df %>% slice(row_indices12)
    
    df6 <- df5 %>%
      slice(-row_indices12:-length(row_names)) %>%
      bind_rows(row11)  %>%  
      bind_rows(row12) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices11, -row_indices12) %>% slice(-1:-(row_indices12-1)))
    
    df <- df6
  }
  
  if(length(levels) == 7){
    raw_row_names <- df$Pathway1
    raw_row <- levels[1:7]
    raw_row_indices <- which(raw_row_names %in% raw_row)
    
    df <- df %>% slice(raw_row_indices)
    
    row_names <- df$Pathway2
    
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
      bind_rows(row3)  %>%  
      bind_rows(row4) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices3, -row_indices4) %>% slice(-1:-(row_indices4-1)))
    
    row_indices5 <- which(row_names %in% levels[3])
    row_indices6 <- which(row_names %in% maxPathway[3])
    
    row5 <- df %>% slice(row_indices5)
    row6 <- df %>% slice(row_indices6)
    
    df3 <- df2 %>%
      slice(-row_indices6:-length(row_names)) %>%
      bind_rows(row5)  %>%  
      bind_rows(row6) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices5, -row_indices6) %>% slice(-1:-(row_indices6-1)))
    
    row_indices7 <- which(row_names %in% levels[4])
    row_indices8 <- which(row_names %in% maxPathway[4])
    
    row7 <- df %>% slice(row_indices7)
    row8 <- df %>% slice(row_indices8)
    
    df4 <- df3 %>%
      slice(-row_indices8:-length(row_names)) %>%
      bind_rows(row7)  %>%  
      bind_rows(row8) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices7, -row_indices8) %>% slice(-1:-(row_indices8-1)))
    
    row_indices9 <- which(row_names %in% levels[5])
    row_indices10 <- which(row_names %in% maxPathway[5])
    
    row9 <- df %>% slice(row_indices9)
    row10 <- df %>% slice(row_indices10)
    
    df5 <- df4 %>%
      slice(-row_indices10:-length(row_names)) %>%
      bind_rows(row9)  %>%  
      bind_rows(row10) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices9, -row_indices10) %>% slice(-1:-(row_indices10-1)))
    
    row_indices11 <- which(row_names %in% levels[6])
    row_indices12 <- which(row_names %in% maxPathway[6])
    
    row11 <- df %>% slice(row_indices11)
    row12 <- df %>% slice(row_indices12)
    
    df6 <- df5 %>%
      slice(-row_indices12:-length(row_names)) %>%
      bind_rows(row11)  %>%  
      bind_rows(row12) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices11, -row_indices12) %>% slice(-1:-(row_indices12-1)))
    
    row_indices13 <- which(row_names %in% levels[7])
    row_indices14 <- which(row_names %in% maxPathway[7])
    
    row13 <- df %>% slice(row_indices13)
    row14 <- df %>% slice(row_indices14)
    
    df7 <- df6 %>%
      slice(-row_indices14:-length(row_names)) %>%
      bind_rows(row13)  %>%  
      bind_rows(row14) %>%
      bind_rows(df %>% slice(1:length(row_names)) %>% slice(-row_indices13, -row_indices14) %>% slice(-1:-(row_indices14-1)))
    
    df <- df7
  }
  
  
  df$Pathway2 <- fct_inorder(df$Pathway2) %>% fct_rev()
  
  barplot <- ggplot(df, aes(x = NES, y = Pathway2, fill = Pathway1))+
    geom_col()+ 
    geom_text(aes(label = Label), 
              color = 'black',
              size = num_size,
              hjust = "left",
              nudge_x = 0.1) + 
    scale_x_continuous(limits = xlim,
                       expand = expansion(mult = c(0, .1))) + 
    labs(x = "NES", y = "", title = title) +
    scale_fill_manual(values = pathway_color)
  
  g <- ggplot_build(barplot)
  mycol <- g$data[[1]]["fill"]
  col <- rev(mycol$fill)
  
  num <- rev(df$NES)
  index <- which(num == 0)
  col[index] <- category_color
  
  
  theme <-  theme_bw()+
    theme(plot.title = element_text(size = rel(1),
                                    hjust = 0.5,
                                    face = 'plain'),
          axis.title = element_text(size = rel(1)),
          axis.text.x = element_text(size = text_x_size,
                                     color = 'black'),
          axis.text.y = element_text(size = text_y_size,
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
