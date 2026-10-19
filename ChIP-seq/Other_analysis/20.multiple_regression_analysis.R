setwd("")
library(dplyr)
library(data.table)
library(ggplot2)
library(pheatmap)
library(broom)

#Multiple regression analysis

#This code is used for the following figure
#Fig. 4E, 4F

#Processing data
data1 <- fread("Table_S3_NGS_data_in_hotspots_1.csv")
data2 <- fread("Table_S3_NGS_data_in_hotspots_3.csv")

pseudocount <- 0.0001
data <- merge(data1, data2, by = "name", sort = FALSE) %>% 
  mutate(Rec12_WT_plus_ave  = (Rec12_WT_plus_rep1  + Rec12_WT_plus_rep2)  / 2,
         Rec12_WT_minus_ave = (Rec12_WT_minus_rep1 + Rec12_WT_minus_rep2) / 2) %>% 
  mutate(
    Rec12_FC   = log2((Rec12_WT_minus_ave    + pseudocount) / (Rec12_WT_plus_ave    + pseudocount)),
    Rec8_FC    = log2((Rec8_WT_minus_rep1    + pseudocount) / (Rec8_WT_plus_rep1    + pseudocount)),
    H3_FC      = log2((H3_WT_minus_rep1      + pseudocount) / (H3_WT_plus_rep1      + pseudocount))
  ) %>% 
  select(Rec12_FC, Rec8_FC, H3_FC)

#Normalize by z-score
data_normalize <- scale(data) %>% 
  as.data.table()

#Multiple regression
fit <- lm(Rec12_FC ~ Rec8_FC * H3_FC,
          data = data_normalize)

#Multiple regression result (Fig. 4E)
result_table <- tidy(fit) %>% 
  as.data.table()


#Prediction of Rec12 FC
FC_range <- seq(-1, 1, length.out = 100)

#x, scaling of Rec8_FC
FC_range_scaled_rec8 <- (FC_range - mean(data$Rec8_FC)) / sd(data$Rec8_FC)
#y, scaling of H3_FC
FC_range_scaled_h3 <- (FC_range - mean(data$H3_FC)) / sd(data$H3_FC)

fit_matrix <- lapply(FC_range_scaled_rec8, function(x){
  sapply(rev(FC_range_scaled_h3), function(y){
    result_table$estimate[1] + result_table$estimate[2] * x + result_table$estimate[3] * y + result_table$estimate[4] * x * y
  })
}) %>% 
  do.call(cbind, .) %>% 
  as.data.frame()

#back from z-score
rec12_mean <- mean(data$Rec12_FC)
rec12_sd   <- sd(data$Rec12_FC)
fit_matrix_rm_scaling <- fit_matrix * rec12_sd + rec12_mean

max_val <- max(abs(fit_matrix_rm_scaling))

#Heat map showing the predicted Rec12 FC (Fig. 4F)
pheatmap(fit_matrix_rm_scaling,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         border_color = NA,
         show_rownames = FALSE,
         show_colnames = FALSE,
         cellwidth  = 3,
         cellheight = 3,
         breaks = seq(-1.5, 1.5, length.out = 101), 
         legend = FALSE)

