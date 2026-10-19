setwd("")
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)

#Violin plot

#This code is used for the following figure
#Fig. 3C, 3E, 3G, S3D, S3F
#Fig. 4C

#function for violin plot
Violin_plot<- function(data,
                       xcol,
                       ycol,
                       xcol_order = NULL,
                       log2 = FALSE, 
                       pseudocount = 0,
                       ylim = NULL, 
                       x_label = NULL, 
                       y_label = NULL, 
                       box_width = 2, 
                       line_width = 1,
                       font_size = 12,
                       plot_margin = c(0.5,0.5,0.5,0.5)){
  #data for plot
  ycol_val <- data[[ycol]]
  
  #log2
  if(isTRUE(log2)){
    y_val <- log2((ycol_val + pseudocount))
  }else{
    y_val <- (ycol_val + pseudocount)
  }
  
  x_val <- data[[xcol]]
  
  #df for plot
  plot_data <- data.frame(
    x_val = x_val,
    y_val = y_val
  )
  if(is.null(xcol_order) == TRUE){
    plot_data$x_val <- factor(plot_data$x_val, levels = unique(plot_data$x_val))
  }else{
    plot_data$x_val <- factor(plot_data$x_val, levels = xcol_order) 
  }
  
  #median
  med <- plot_data %>%
    group_by(x_val) %>%
    summarise(median = median(y_val, na.rm = TRUE))
  print(med)
  
  final_fig <- ggplot()+
    geom_violin(data  = plot_data, aes(x = x_val, y = y_val,  fill = x_val))+ 
    geom_boxplot(data = plot_data, aes(x = x_val, y = y_val), fill = "white", width = box_width, linewidth = line_width, outlier.shape = NA)+
    # geom_beeswarm(data = plot_data, aes(x = cluctering, y = y_val, fill = cluctering), size = 1, cex = 2)+
    scale_y_continuous(limits = ylim, expand = c(0,0), labels = scales::label_number(accuracy = 0.1))+
    labs(x = x_label,
         y = y_label, 
         fill  = NULL, 
         title = NULL)+
    theme_classic(base_size = font_size)+
    theme(
      text        = element_text(size = font_size, color = "black"),
      axis.text   = element_text(size = font_size, color = "black"),
      axis.text.x = element_blank(),
      legend.text = element_text(size = font_size, color = "black"),
      plot.title  = element_text(size = font_size, color = "black", hjust = 0.5),
      legend.position = "none",
      plot.margin = unit(plot_margin, "cm")
    )
  
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
xcol_name_ls  <- "clustering"
xcol_order_ls <- c("increased", "ns", "decreased")
ycol_name_ls  <- c("Rec8_FC", "H3_FC", "Atf1_FC", "H3K9ac_FC", "H3K4me3_FC", "RNA_FC")

ylim_ls <- list(
  c(-1.5, 1.5), c(-2.5, 2), c(-1.5, 2.5), c(-2, 4), c(-2, 3), c(-4,12)
)


#Figure
#log2 and pseudocount are already considered during data processing
Fig_ls <- lapply(1:length(ycol_name_ls), function(n){
  
  Violin_plot(data = data, 
              xcol = xcol_name_ls, 
              ycol = ycol_name_ls[n], 
              xcol_order = xcol_order_ls,
              log2 = FALSE, 
              pseudocount = 0,
              ylim = ylim_ls[[n]], 
              x_label = NULL, 
              y_label = paste0(ycol_name_ls[n], " (log2)"), 
              box_width = 0.3, 
              line_width = 0.5, 
              font_size = 10, 
              plot_margin = c(0.5,0.5,0.5,0.5))
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
