setwd("")
library(data.table)
library(dplyr)
library(limma)

#statistical analysis of Rec12 ChIP-seq using limma
#resulting p-value are indicated in "Table_S5_NGS_data_in_hotspots.xlsx"

#e.g., Glucose+ vs Glucose- in WT (Fig. 2B)
data <- fread("Table_S5_NGS_data_in_hotspots_1.csv")

sample_col <- c("Rec12_WT_plus_rep1",
                "Rec12_WT_plus_rep2",
                "Rec12_WT_minus_rep1",
                "Rec12_WT_minus_rep2")

pseudocount <- 0.0001
df_for_limma <- data %>% 
  select(all_of(sample_col)) 
df_for_limma <- log2(df_for_limma + pseudocount)

group  <- factor(c("plus", "plus", "minus", "minus"), levels = c("plus", "minus"))
design <- model.matrix(~ group)
fit    <- lmFit(df_for_limma, design)
fit2   <- eBayes(fit, trend = FALSE, robust = FALSE)

limma_result <- topTable(fit2, coef = "groupminus", number = Inf, adjust.method = "fdr", sort.by = "none") %>%
  mutate(significance = case_when(
    adj.P.Val <  0.05 ~ "significant",
    adj.P.Val >= 0.05 ~ "non-significant"
  )) %>% 
  select(adj.P.Val, significance)


#e.g., FC (Glucose-/Glucose+) in WT vs FC in pcr1d
data <- fread("Table_S5_NGS_data_in_hotspots_2.csv")
valid_index  <- which(data$mask == "not_masked")
data2 <- data %>% 
  filter(mask == "not_masked")

sample_col <- c("Rec12_WT_plus_rep1",
                "Rec12_WT_plus_rep2",
                "Rec12_WT_plus_rep3",
                "Rec12_WT_minus_rep1",
                "Rec12_WT_minus_rep2",
                "Rec12_WT_minus_rep3",
                "Rec12_pcr1d_plus_9h_rep1",
                "Rec12_pcr1d_plus_9h_rep2",
                "Rec12_pcr1d_minus_9h_rep1",
                "Rec12_pcr1d_minus_9h_rep2")

pseudocount <- 0.0001
df_for_limma <- data2 %>% 
  select(all_of(sample_col)) 
df_for_limma <- log2(df_for_limma + pseudocount)

group  <- factor(c(rep(c("plus", "minus"), each = 3),
                   rep(c("plus", "minus"), each = 2)), 
                 levels = c("plus", "minus"))

genotype <- factor(c(rep("WT", times = 6), rep("pcr1d", times = 4)),
                   levels = c("WT", "pcr1d"))

design <- model.matrix(~ genotype * group)

fit    <- lmFit(df_for_limma, design)
fit2   <- eBayes(fit, trend = FALSE, robust = FALSE)

limma_result <- topTable(fit2, coef = "genotypepcr1d:groupminus", number = Inf, adjust.method = "fdr", sort.by = "none") %>%
  mutate(significance = case_when(
    adj.P.Val <  0.05 ~ "significant",
    adj.P.Val >= 0.05 ~ "non-significant"
  )) %>% 
  select(adj.P.Val, significance)

limma_result_2 <- data.table(
  adj.P.Val    = rep(NA, times = nrow(data)),
  significance = rep(NA, times = nrow(data))
)
limma_result_2$adj.P.Val[valid_index]    <- limma_result$adj.P.Val
limma_result_2$significance[valid_index] <- limma_result$significance


