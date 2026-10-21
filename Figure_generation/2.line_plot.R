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
                      summary_group_col = "condition",
                      line_group_col    = "clustering",
                      xcol = "time",
                      ycol = "value", 
                      xlab = "", 
                      ylab = "",
                      font_size = 18){
  
  current_data <- data
  
  names(current_data)[names(current_data) == summary_group_col] <- "sample_info"
  condition_name <- unique(current_data$sample_info)
  current_data$sample_info <- factor(current_data$sample_info, 
                                     levels = condition_name)
  
  names(current_data)[names(current_data) == line_group_col] <- "clustering"
  clustering_name <- unique(current_data$clustering)
  current_data$clustering <- factor(current_data$clustering, 
                                    levels = clustering_name)
  
  names(current_data)[names(current_data) == xcol] <- "xcol"
  names(current_data)[names(current_data) == ycol] <- "ycol"
  
  print(current_data)
  
  average_df <- current_data %>%
    group_by(sample_info) %>%
    summarise(
      clustering = first(clustering),
      xcol       = first(xcol),
      average    = mean(ycol, na.rm = TRUE),
      SD         = sd(ycol, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    arrange(xcol)
  
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
  mutate(summary_group = paste0(strain, "_", condition, "_", time, "h"),
         line_group    = paste0(strain, "_", condition))

Fig_6A <- line_plot(data = data2,
                    xlim = c(3.5,12.5),
                    ylim = c(0,100),
                    summary_group_col = "summary_group", 
                    xcol = "time",
                    ycol = "N2-4_ratio",
                    line_group_col = "line_group",
                    xlab = "Time",
                    ylab = "MI + MII")
Fig_6A


