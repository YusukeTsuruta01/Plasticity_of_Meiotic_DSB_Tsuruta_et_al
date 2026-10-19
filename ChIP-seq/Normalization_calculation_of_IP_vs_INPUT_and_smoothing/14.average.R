setwd("")
library(data.table)
library(dplyr)
library(rtracklayer)

#data
base_path <- ""
path <- paste0(base_path, "/coverage_data/smoothed_data/501bp/RData/")
output_path <- paste0(base_path, "/coverage_data/average/")

#501bp smoothed
window_size <- 501
load(paste0(path, 
            "/Rec12_IP_vs_INPUT_",
            window_size, 
            "_smoothed.RData.gz"))

rep1 <- c(1,2)
rep2 <- c(3,4)

output_name <- c(
  ""
)


Rec12_IP_vs_INPUT_smoothed_ave <- lapply(1:length(rep1), function(n){
  current_data_rep1 <- Rec12_IP_vs_INPUT_smoothed[[rep1[n]]]
  current_data_rep2 <- Rec12_IP_vs_INPUT_smoothed[[rep2[n]]]
  
  result_ls <- lapply(names(current_data_rep1), function(i){
    result <- current_data_rep1[[i]]
    ave    <- (current_data_rep1[[i]]$count + current_data_rep2[[i]]$count) / 2
    result$count <- ave
    return(result)
  })
  names(result_ls) <- names(current_data_rep1)
  return(result_ls)
})
names(Rec12_IP_vs_INPUT_smoothed_ave) <- output_name

dir.create(output_path)
dir.create(paste0(output_path, "/", window_size, "bp/"))
dir.create(paste0(output_path, "/", window_size, "bp/RData/"))

save(Rec12_IP_vs_INPUT_smoothed_ave, 
     file = paste0(output_path, 
                   "/",
                   window_size, 
                   "bp/RData/Rec12_IP_vs_INPUT_", 
                   window_size, 
                   "_smoothed_average.RData.gz"), compress = TRUE)

#for bigwig
change_format <- function(data_ls, used_col){
  data_df <- do.call(rbind, data_ls)
  data_df[[used_col]][is.na(data_df[[used_col]])] <- 0
  
  granges <- GRanges(
    seqnames = data_df$chr,
    ranges   = IRanges(start = data_df$coordinate,
                       end   = data_df$coordinate),
    score    = data_df[[used_col]]
  )
  
  seqinfo(granges) <- Seqinfo(
    seqnames   = c("I",     "II",    "III",  "AB325691", "MT", "MTR"),
    seqlengths = c(5579133, 4539804, 2452883, 20000,     19431, 20128)
  )
  return(granges)
}
Rec12_IP_vs_INPUT_smoothed_ave_bigwig <- lapply(Rec12_IP_vs_INPUT_smoothed_ave, function(n){
  change_format(data_ls  = n, 
                used_col = "count"
  )
})

#save as bigwig
dir.create(paste0(output_path, "/", window_size, "bp/bigwig/"))

for (i in 1:length(Rec12_IP_vs_INPUT_smoothed_ave_bigwig)) {
  export(Rec12_IP_vs_INPUT_smoothed_ave_bigwig[[i]], 
         paste0(output_path, 
                "/", 
                window_size, 
                "bp/bigwig/Rec12_IP_vs_INPUT_", 
                names(Rec12_IP_vs_INPUT_smoothed_ave_bigwig)[i], 
                "_", 
                window_size, 
                "bp_smoothed_average.bw"), 
         format = "BigWig")
}
