setwd("")
library(dplyr)
library(data.table)
library(FSA)

#statistical analysis
#Kruskal Wallis and Dunn

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

test_col_name_ls  <- c("Rec8_FC", "H3_FC", "Atf1_FC", "H3K9ac_FC", "H3K4me3_FC", "RNA_FC")
data_ls <- lapply(test_col_name_ls, function(n){
  result <- data %>% 
    select(!!sym(n), clustering)
  names(result)[1] <- "value"
  return(result)
})
names(data_ls) <- test_col_name_ls

#Kruskal Wallis
kruskal_ls <- lapply(data_ls, function(n){
  kruskal.test(value ~ clustering, data = n)
})

kruskal_p_ls <- sapply(kruskal_ls, function(n){
  n$p.value
})

#Dunn
dunn_ls <- lapply(data_ls, function(n){
  dunnTest(value ~ clustering, data = n, method = "bh")
})

dunn_p_ls <- lapply(dunn_ls, function(n){
  n$res
})

dunn_p_total <- do.call(rbind, dunn_p_ls) %>% 
  as.data.table() %>% 
  mutate(term      = rep(test_col_name_ls, each = 3),
         Kruskal_p = rep(kruskal_p_ls,     each = 3))



