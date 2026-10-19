setwd("")
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)

#Scatter plot

#This code is used for the following figure
#Fig. S1A, S1E
#Fig. 3B, 3D, 3F, S3B, S3C, S3E
#Fig. 4B, 4D, S4C, S4D, 
#Fig. 5C, 5D, 5F, S5A, S5B, S5C, S5G
#Fig. 6D, 6F, S6C, S6D

#function for scatter plot
Scatter_plot <- function(data,
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
                         y_x_line = FALSE){
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
  
  #correlation
  correlation  <- cor.test(plot_data$x_val,
                           plot_data$y_val,
                           use = "complete.obs",
                           method = "pearson")
  r_value_text   <- paste0("r = ", signif(as.numeric(correlation$estimate), digits = 3))
  r_value_p_text <- paste0("p = ", signif(as.numeric(correlation$p.value),  digits = 3))
  
  #linear regression
  linear_regression <- lm(y_val ~ x_val, data = plot_data)
  #print(linear_regression)
  linear_regression_1 <- paste0("y = ", 
                                signif(as.numeric(coef(linear_regression)[2]), digits = 3), 
                                "x",
                                ifelse(coef(linear_regression)[1] > 0, "+", ""),
                                signif(as.numeric(coef(linear_regression)[1]), digits = 3))
  
  
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
    geom_smooth(data = plot_data, aes(x = x_val, y = y_val), method = "lm", se = FALSE, colour = "black", linewidth = 0.5)+
    annotate("text", x = max(xlim)*0.2, y = max(ylim)*0.9, label = r_value_text,        size = font_size/3, color = "black")+
    annotate("text", x = max(xlim)*0.2, y = max(ylim)*0.8, label = r_value_p_text,      size = font_size/3, color = "black")+
    annotate("text", x = max(xlim)*0.2, y = max(ylim)*0.7, label = linear_regression_1, size = font_size/3, color = "black")+
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
      legend.position = "none",
      plot.margin = unit(plot_margin, "cm")
    )
  
  if(isTRUE(fix_fig_ratio)){
    final_fig <- final_fig+
      coord_fixed(ratio = 1)
  }
  
  if(isTRUE(y_x_line)){
    final_fig <- final_fig+
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5)
  }
  
  return(final_fig)
}

#e.g.,
#Fig. 3
#Processing data
data1 <- fread("Table_S3_NGS_data_in_hotspots_1.csv")
data2 <- fread("Table_S3_NGS_data_in_hotspots_3.csv")
data3 <- fread("Table_S3_NGS_data_in_hotspots_4.csv")

pseudocount <- 0.0001
data <- merge(data1, data2, by = "name", sort = FALSE) %>% 
  merge(., data3, by = "name", sort = FALSE) %>% 
  mutate(Rec12_WT_plus_ave  = (Rec12_WT_plus_rep1  + Rec12_WT_plus_rep2)  / 2,
         Rec12_WT_minus_ave = (Rec12_WT_minus_rep1 + Rec12_WT_minus_rep2) / 2) %>% 
  mutate(
    Rec12_FC   = log2((Rec12_WT_minus_ave    + pseudocount) / (Rec12_WT_plus_ave    + pseudocount)),
    Rec8_FC    = log2((Rec8_WT_minus_rep1    + pseudocount) / (Rec8_WT_plus_rep1    + pseudocount)),
    H3_FC      = log2((H3_WT_minus_rep1      + pseudocount) / (H3_WT_plus_rep1      + pseudocount)),
    Atf1_FC    = log2((Atf1_WT_minus_rep1    + pseudocount) / (Atf1_WT_plus_rep1    + pseudocount)),
    H3K9ac_FC  = log2((H3K9ac_ratio_WT_minus_rep1  + pseudocount) / (H3K9ac_ratio_WT_plus_rep1  + pseudocount)),
    H3K4me3_FC = log2((H3K4me3_ratio_WT_minus_rep1 + pseudocount) / (H3K4me3_ratio_WT_plus_rep1 + pseudocount)),
    RNA_FC     = log2((TPM_WT_minus_rep1     + pseudocount) / (TPM_WT_plus_rep1     + pseudocount))
  ) %>% 
  mutate(clustering = case_when(
    adj.P.Val < 0.05  & Rec12_FC > 0 ~ "increased",
    adj.P.Val < 0.05  & Rec12_FC < 0 ~ "decreased",
    adj.P.Val >= 0.05                ~ "ns"
  ))
# table(data$clustering)


#parameters for plot
xcol_name_ls <- c("Rec8_FC", "H3_FC", "Atf1_FC", "H3K9ac_FC", "H3K4me3_FC", "RNA_FC")
ycol_name_ls <- "Rec12_FC"

xlim_ls <- list(
  c(-1.5, 1.5), c(-2.5, 2), c(-1.5, 2.5), c(-2, 4), c(-2, 3), c(-4,12)
)
ylim_ls <- c(-3.5,3)


#Figure
#log2 and pseudocount are already considered during data processing
Fig_ls <- lapply(1:length(xcol_name_ls), function(n){
  
  Scatter_plot(data = data, 
               xcol = xcol_name_ls[n], 
               ycol = ycol_name_ls, 
               log2 = FALSE, 
               pseudocount = 0,
               clustering = TRUE,
               clustering_col = "clustering", 
               xlim = xlim_ls[[n]], 
               ylim = ylim_ls, 
               x_label = paste0(xcol_name_ls[n], "(log2)"), 
               y_label = paste0(ycol_name_ls,    "(log2)"), 
               dot_size = 0.5, 
               alpha = 1, 
               font_size = 10, 
               plot_margin = c(0.5,0.5,0.5,0.5), 
               fix_fig_ratio = FALSE, 
               y_x_line = FALSE)
})

all_fig <- cowplot::align_plots(
  Fig_ls[[1]], Fig_ls[[2]], Fig_ls[[3]], Fig_ls[[4]], Fig_ls[[5]], Fig_ls[[6]],
  align = "hv",
  axis  = "tblr"
)

Fig <- plot_grid(
  plotlist = all_fig,
  nrow = 2
)
Fig


