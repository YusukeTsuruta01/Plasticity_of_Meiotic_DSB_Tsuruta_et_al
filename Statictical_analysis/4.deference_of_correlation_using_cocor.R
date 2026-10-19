setwd("")
library(dplyr)
library(data.table)
library(cocor)

#statistical analysis to analyze the difference of correlation using cocor.
#e.g.,
#Fig. 5F
#Processing data
data1 <- fread("Table_S3_NGS_data_in_hotspots_6.csv")
data2 <- fread("Table_S3_NGS_data_in_hotspots_7.csv")

pseudocount <- 0.0001
data <- merge(data1, data2, by = "name", sort = FALSE) %>% 
  filter(mask == "not_masked") %>% 
  mutate(Rec12_WT_plus_ave     = (Rec12_WT_plus_rep1  + Rec12_WT_plus_rep2  + Rec12_WT_plus_rep3)  / 3,
         Rec12_WT_minus_ave    = (Rec12_WT_minus_rep1 + Rec12_WT_minus_rep2 + Rec12_WT_minus_rep3) / 3,
         Rec12_rec8d_plus_ave  = (Rec12_rec8d_plus_rep1  + Rec12_rec8d_plus_rep2)  / 2,
         Rec12_rec8d_minus_ave = (Rec12_rec8d_minus_rep1 + Rec12_rec8d_minus_rep2) / 2
         ) %>% 
  mutate(
    Rec12_WT_FC    = log2((Rec12_WT_minus_ave    + pseudocount) / (Rec12_WT_plus_ave    + pseudocount)),
    Rec12_rec8d_FC = log2((Rec12_rec8d_minus_ave + pseudocount) / (Rec12_rec8d_plus_ave + pseudocount)),
    H3_FC          = log2((H3_WT_minus_rep1      + pseudocount) / (H3_WT_plus_rep1      + pseudocount))
  )
# %>% 
  # select(name, Rec12_WT_FC, Rec12_rec8d_FC, H3_FC)

#cocor
cocor.dep.groups.overlap(r.jk = cor(data$Rec12_WT_FC,    data$H3_FC,          use = "complete.obs", method = "pearson"),
                         r.jh = cor(data$Rec12_rec8d_FC, data$H3_FC,          use = "complete.obs", method = "pearson"),
                         r.kh = cor(data$Rec12_WT_FC,    data$Rec12_rec8d_FC, use = "complete.obs", method = "pearson"), 
                         n = nrow(data))

# meng1992: Meng, Rosenthal, and Rubin's z (1992)
#   z = 2.4872, p-value = 0.0129
#   Null hypothesis rejected
#   95% confidence interval for r.jk - r.jh: 0.0227 0.1913
#   Null hypothesis rejected (Interval does not include 0)
