setwd("")
library(dplyr)
library(data.table)

#statistical analysis
#ANOVA and Tukey HSD
#resulting p-value are indicated in "Table_S5_NGS_data_in_hotspots.xlsx"

#e.g., Fig.6E
data <- fread("Table_S1_raw_data_Fig.6E.csv")

#use only 2 h data
data2 <- data %>% 
  filter(time == 2) %>% 
  mutate(condition = paste0(strain, "_", sub("Glucose_","", condition), "_", time, "h"))

gene_ls <- c("osh6", "SPBC2A9.02")
data_ls <- lapply(gene_ls, function(n){
  data2 %>% 
    filter(gene == n)
})

#ANOVA (one-way)
anova_ls <- lapply(data_ls, function(n){
  aov(ratio ~ condition, data = n)
})
anova_p_ls <- lapply(anova_ls, function(n){
  summary(n)
})

#TukeyHSD
Tukey_result <- lapply(anova_ls, function(n){
  TukeyHSD(n) 
})




