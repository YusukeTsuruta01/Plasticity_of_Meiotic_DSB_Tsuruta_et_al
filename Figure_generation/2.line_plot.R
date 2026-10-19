setwd("")
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggbeeswarm)

#Line plot

#This code is used for the following figure
#Fig. S2E
#Fig. 6A

#function for line plot
line_plot <- function(data,
                      xlim,
                      ylim, 
                      condition_col  = "condition",
                      clustering_col = "clustering",
                      xcol = "time",
                      ycol = "value", 
                      xlab = "", 
                      ylab = "",
                      font_size = 18){
  
  current_data <- data
  
  names(current_data)[names(current_data) == condition_col] <- "sample_info"
  condition_name <- unique(current_data$sample_info)
  current_data$sample_info <- factor(current_data$sample_info, 
                                     levels = condition_name)
  
  names(current_data)[names(current_data) == clustering_col] <- "clustering"
  clustering_name <- unique(current_data$clustering)
  current_data$clustering <- factor(current_data$clustering, 
                                    levels = clustering_name)
  
  names(current_data)[names(current_data) == xcol] <- "xcol"
  names(current_data)[names(current_data) == ycol] <- "ycol"
  
  # print(current_data)
  
  average <- sapply(condition_name, function(n){
    current_data %>% 
      filter(sample_info == n) %>% 
      pull(ycol) %>% 
      mean(na.rm = TRUE)
  })
  
  SD <- sapply(condition_name, function(n){
    current_data %>% 
      filter(sample_info == n) %>% 
      pull(ycol) %>% 
      sd(na.rm = TRUE)
  })
  
  xcol <- sapply(condition_name, function(n){
    current_data %>% 
      filter(sample_info == n) %>% 
      pull(xcol) %>% 
      unique()
  })
  
  clustering_name_for_ave_df <- sapply(condition_name, function(n){
    current_data %>% 
      filter(sample_info == n) %>% 
      pull(clustering) %>% 
      unique()
  })
  
  average_df <- data.frame(
    clustering  = clustering_name_for_ave_df,
    xcol        = xcol,
    average     = average,
    SD          = SD
  )
  
  average_df$clustering <- factor(average_df$clustering, 
                                  levels = unique(average_df$clustering))
  
  # print(average_df)
  
  ggplot()+
    geom_line(data  = average_df,   aes(x = xcol, y = average, colour = clustering), linewidth = 1)+
    geom_point(data = current_data, aes(x = xcol, y = ycol,    colour = clustering), size = 2)+
    # geom_errorbar(data = average_df, aes(x = xcol, ymin = ycol - SD, ymax = ycol + SD), width = 0.2)+
    scale_x_continuous(limits = xlim, expand = c(0,0), labels = scales::label_number(accuracy = 0.1))+
    scale_y_continuous(limits = ylim, expand = c(0,0), labels = scales::label_number(accuracy = 0.1))+
    labs(x = xlab, 
         y = ylab, 
         fill = NULL,
         title = NULL)+
    theme_classic(base_size = font_size)+
    theme(
      text        = element_text(size = font_size, color = "black"),
      axis.text   = element_text(size = font_size, color = "black"),
      plot.title  = element_text(size = font_size, color = "black", hjust = 0.5)
    )
}

#e.g.,
#Fig. 6A
data <- fread("Table_S1_raw_data_Fig.6A.csv")
data2 <- data %>% 
  mutate(condition2 = paste0(strain, "_", condition, "_", time, "h"),
         condition3 = paste0(strain, "_", condition))

Fig_6A <- line_plot(data = data2,
                    xlim = c(3.5,12.5),
                    ylim = c(0,100),
                    condition_col = "condition2", 
                    xcol = "time",
                    ycol = "N2-4_ratio",
                    clustering_col = "condition3",
                    xlab = "Time",
                    ylab = "MI + MII")
Fig_6A


