setwd("")
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(pheatmap)

#Heat map of correlation matrix

#This code is used for the following figure
#Fig. 1E, S1D
#Fig. S3G
#Fig. S4F

#e.g.,
#Fig. S3G

#Processing data
data1 <- fread("Table_S3_NGS_data_in_hotspots_1.csv")
data2 <- fread("Table_S3_NGS_data_in_hotspots_3.csv")
data3 <- fread("Table_S3_NGS_data_in_hotspots_4.csv")

data <- merge(data1, data2, by = "name", sort = FALSE) %>% 
  merge(., data3, by = "name", sort = FALSE) %>% 
  mutate(Rec12_WT_plus_ave  = (Rec12_WT_plus_rep1  + Rec12_WT_plus_rep2)  / 2,
         Rec12_WT_minus_ave = (Rec12_WT_minus_rep1 + Rec12_WT_minus_rep2) / 2)

sample_name <- c("Rec12", "Atf1", "H3", "H3K9ac", "H3K4me3", "RNA")
plus_name   <- c("Rec12_WT_plus_ave",  "Atf1_WT_plus_rep1",  "H3_WT_plus_rep1",  "H3K9ac_ratio_WT_plus_rep1",  "H3K4me3_ratio_WT_plus_rep1",  "TPM_WT_plus_rep1")
minus_name  <- c("Rec12_WT_minus_ave", "Atf1_WT_minus_rep1", "H3_WT_minus_rep1", "H3K9ac_ratio_WT_minus_rep1", "H3K4me3_ratio_WT_minus_rep1", "TPM_WT_minus_rep1")

#calculation of FC
FC_df <- lapply(1:length(plus_name), function(n){
    log2((data[[minus_name[n]]] + 0.0001) / (data[[plus_name[n]]] + 0.0001))
  }) %>% 
  do.call(cbind, .) %>% 
  as.data.table()
names(FC_df) <- sample_name

#calculation of correlation
correlation_matrix <- lapply(1:ncol(FC_df), function(n){
  sapply(1:ncol(FC_df), function(i){
    cor(x = FC_df[[i]],
        y = FC_df[[n]],
        use = "complete.obs")
    })
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

#plot
pheatmap(correlation_matrix, 
         cluster_cols = FALSE, 
         cluster_rows = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-1, 1, length.out = 101),
         display_numbers = TRUE,
         number_format = "%.3g",
         main = "correlation of ChIP-seq",
         fontsize = 12,
         fontsize_number = 12,
         labels_col = total_names_ls,
         labels_row = total_names_ls,
         angle_col = 45,
         cellwidth = 60,
         cellheight = 60
)
