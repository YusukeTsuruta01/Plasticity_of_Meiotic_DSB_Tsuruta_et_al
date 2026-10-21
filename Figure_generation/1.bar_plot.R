setwd("")
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggbeeswarm)

#Bar plot

#This code is used for the following figure
#Fig. 1B
#Fig. 2C, 2E
#Fig. S2C, S2D
#Fig. 5A
#Fig. 6B, 6E, 6G
#Fig. S6A, S6B

#function for bar plot
bar_plot <- function(data, 
                     ylim, 
                     summary_group_col = "condition",
                     result_col = "value", 
                     ylab = "", 
                     font_size = 18){
  
  current_data <- data
  
  names(current_data)[names(current_data) == summary_group_col] <- "sample_info"
  condition_name <- unique(current_data$sample_info)
  current_data$sample_info <- factor(current_data$sample_info, 
                                     levels = condition_name)
  
  names(current_data)[names(current_data) == result_col] <- "value"
  
  # print(current_data)
  
  average_df <- current_data %>%
    group_by(sample_info) %>%
    summarise(
      average = mean(value, na.rm = TRUE),
      SD      = sd(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  average_df$sample_info <- factor(average_df$sample_info, 
                                   levels = unique(average_df$sample_info))
  # print(average_df)
  
  ggplot()+
    geom_col(data = average_df, aes(x = sample_info, y = average, fill = sample_info), width = 0.6)+
    geom_errorbar(data = average_df, aes(x = sample_info, ymin = average - SD, ymax = average + SD), width = 0.2)+
    geom_beeswarm(data = current_data, aes(x = sample_info, y = value), fill = "white", shape = 21, colour = "black", size = 3)+
    scale_y_continuous(limits = ylim, expand = c(0,0), labels = scales::label_number(accuracy = 0.1))+
    labs(x = NULL, 
         y = ylab, 
         fill = NULL,
         title = NULL)+
    theme_classic(base_size = font_size)+
    theme(
      text        = element_text(size = font_size, color = "black"),
      axis.text.x = element_blank(),
      # axis.text.x = element_text(size = font_size, color = "black", angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = font_size, color = "black"), 
      legend.text = element_text(size = font_size, color = "black"), 
      plot.title  = element_text(size = font_size, color = "black", hjust = 0.5)
    )
}

#e.g.,
#Fig. 6B
data <- fread("Table_S1_raw_data_Fig.6B.csv")
data2 <- data %>% 
  mutate(summary_group = paste0(strain, "_", condition, "_", time, "h"))

Fig_6B <- bar_plot(data = data2,
                   ylim = c(0,8000000),
                   summary_group_col = "summary_group", 
                   result_col        = "total_signal", 
                   ylab = "Total Rec12")
Fig_6B


