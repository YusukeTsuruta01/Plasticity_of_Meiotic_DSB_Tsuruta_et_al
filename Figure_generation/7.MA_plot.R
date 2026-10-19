setwd("")
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)

#MA plot

#This code is used for the following figure
#Fig. 2B, S2B
#Fig. 6D

#function for MA plot
MA_plot <- Scatter_plot <- function(data,
                                    xcol,
                                    ycol,
                                    log2 = FALSE, 
                                    pseudocount = 0,
                                    clustering = FALSE,
                                    clustering_col = NULL, 
                                    xlim = NULL, 
                                    ylim = NULL, 
                                    x_label = NULL, 
                                    y_label = NULL, 
                                    dot_size = 2, 
                                    alpha = 1,
                                    font_size = 12,
                                    plot_margin = c(0.5,0.5,0.5,0.5),
                                    fix_fig_ratio = FALSE, 
                                    med_line = FALSE){
  #data for plot
  xcon_val <- data[[xcol]]
  ycon_val <- data[[ycol]]
  
  #log2
  if(isTRUE(log2)){
    x_val <- log2(xcon_val + pseudocount)
    y_val <- log2(ycon_val + pseudocount)
  }else{
    x_val <- xcon_val + pseudocount
    y_val <- ycon_val + pseudocount
  }
  
  #df for plot
  plot_data <- data.frame(
    x_val = x_val,
    y_val = y_val
  )
  
  #median
  med <- median(plot_data$y_val)
  med_text <- paste0("log2 median = ", signif(med), digits = 3)
  # print(med)
  
  if(isTRUE(clustering)){
    #clustering
    plot_data <- plot_data %>% 
      mutate(cluctering = data[[clustering_col]])
    plot_data$cluctering <- factor(plot_data$cluctering, levels = unique(plot_data$cluctering))
    
    fig <- ggplot()+
      geom_point(data = plot_data, aes(x = x_val, y = y_val, colour = cluctering), size = dot_size, alpha = alpha)
    
  }else{
    fig <- ggplot()+
      geom_point(data = plot_data, aes(x = x_val, y = y_val), size = dot_size)
  }
  
  final_fig <- fig+
    annotate("text", x = max(xlim)*0.2, y = max(ylim)*0.8, label = med_text, size = font_size/3, color = "black")+
    scale_x_continuous(limits = xlim, expand = c(0,0), labels = scales::label_number(accuracy = 0.1), oob = scales::squish)+
    scale_y_continuous(limits = ylim, expand = c(0,0), labels = scales::label_number(accuracy = 0.1), oob = scales::squish)+
    labs(x = x_label, 
         y = y_label, 
         colour = NULL, 
         title  = NULL)+
    theme_classic(base_size = font_size)+
    theme(
      text        = element_text(size = font_size, color = "black"),
      axis.text   = element_text(size = font_size, color = "black"),
      legend.text = element_text(size = font_size, color = "black"),
      plot.title  = element_text(size = font_size, color = "black", hjust = 0.5),
      legend.position = "right",
      plot.margin = unit(plot_margin, "cm")
    )
  
  if(isTRUE(fix_fig_ratio)){
    final_fig <- final_fig+
      coord_fixed(ratio = 1)
  }
  
  if(isTRUE(med_line)){
    final_fig <- final_fig+
      geom_hline(yintercept = med, linetype = "dashed", color = "black", linewidth = 0.5)
  }
  
  return(final_fig)
}

#e.g.,
#Fig. 2B
#Processing data
data <- fread("Table_S3_NGS_data_in_hotspots_1.csv")

pseudocount <- 0.0001
data <- data %>% 
  mutate(Rec12_WT_plus_ave  = (Rec12_WT_plus_rep1  + Rec12_WT_plus_rep2)  / 2,
         Rec12_WT_minus_ave = (Rec12_WT_minus_rep1 + Rec12_WT_minus_rep2) / 2) %>% 
  mutate(
    Rec12_FC  = log2((Rec12_WT_minus_ave + pseudocount) / (Rec12_WT_plus_ave + pseudocount)),
    Rec12_ave = log2((Rec12_WT_plus_ave + Rec12_WT_minus_ave + pseudocount) / 2)
  ) %>% 
  mutate(clustering = case_when(
    adj.P.Val < 0.05  & Rec12_FC > 0 ~ "increased",
    adj.P.Val < 0.05  & Rec12_FC < 0 ~ "decreased",
    adj.P.Val >= 0.05                ~ "ns"
  ))
# table(data$clustering)


#parameters for plot
xcol_name_ls <- "Rec12_ave"
ycol_name_ls <- "Rec12_FC"

xlim_ls <- c(0, 20)
ylim_ls <- c(-3.5,3)


#Figure
#log2 and pseudocount are already considered during data processing
Fig.2B <- Scatter_plot(data = data, 
                       xcol = xcol_name_ls, 
                       ycol = ycol_name_ls, 
                       log2 = FALSE, 
                       pseudocount = 0,
                       clustering = TRUE,
                       clustering_col = "clustering", 
                       xlim = xlim_ls, 
                       ylim = ylim_ls, 
                       x_label = paste0(xcol_name_ls, "(log2)"), 
                       y_label = paste0(ycol_name_ls, "(log2)"), 
                       dot_size = 1, 
                       alpha = 1, 
                       font_size = 10, 
                       plot_margin = c(0.5,0.5,0.5,0.5), 
                       fix_fig_ratio = FALSE, 
                       med_line = TRUE)
Fig.2B

