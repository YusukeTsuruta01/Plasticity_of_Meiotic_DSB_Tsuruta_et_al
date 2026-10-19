setwd("")
library(data.table)
library(dplyr)
library(rtracklayer)

#data
base_path <- ""
path <- paste0(base_path, "/coverage_data/normalized/RData/IP_vs_INPUT/")
output_path <- paste0(base_path, "/coverage_data/smoothed_data/")

smoothing_window_ls <- c(251, 501, 751, 1001)

load(paste0(path, "/Rec12_ChIP_seq_IP_vs_INPUT_10bp_bin_masked.RData.gz"))


#function for smoothing
smoothing <- function(data_ls, window_size){
  parzen.flt <- function(n=51) {
    if (n == 1) 
      c <- 1
    else {
      c <- 1-abs(((0:(n-1))-.5*(n-1))/(.5*(n+1)))
    }
    return(c/sum(c))
  }
  
  result_ls <- lapply(data_ls, function(n){
    current_data <- n
    
    # #na.approx() to complement NA
    # #This is for showing the profile data, do not use this data for any analysis
    # current_data$count <- na.approx(current_data$count, na.rm = FALSE)
    
    smoothed <- current_data$count %>% 
      stats::filter(parzen.flt(n = window_size)) %>% 
      as.numeric()
    current_data$count <- smoothed
    return(current_data)
  })
  names(result_ls) <- names(data_ls)
  return(result_ls)
}

smoothed_ls <- lapply(smoothing_window_ls, function(n){
  lapply(IP_vs_INPUT_normalized_peak_SK1_ls, function(i){
    smoothing(data_ls     = i, 
              window_size = n)
  })
})
names(smoothed_ls) <- paste0("smoothed_", smoothing_window_ls, "bp")

#save as RData
dir.create(output_path)
for (i in 1:length(smoothing_window_ls)) {
  dir.create(paste0(output_path, smoothing_window_ls[i], "bp/"))
  dir.create(paste0(output_path, smoothing_window_ls[i], "bp/RData/"))
  
  Rec12_IP_vs_INPUT_smoothed <- smoothed_ls[[i]]
  save(Rec12_IP_vs_INPUT_smoothed, 
       file = paste0(output_path, 
                     smoothing_window_ls[i], 
                     "bp/RData/Rec12_IP_vs_INPUT_", 
                     smoothing_window_ls[i], 
                     "_smoothed.RData.gz"))
}


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
smoothed_ls_bigwig <- lapply(smoothed_ls, function(n){
  lapply(n, function(i){
    change_format(data_ls  = i, 
                  used_col = "count")
  })
})

#save as bigwig
for (i in 1:length(smoothing_window_ls)) {
  dir.create(paste0(output_path, smoothing_window_ls[i], "bp/bigwig/"))
  
  for (j in 1:length(smoothed_ls_bigwig[[i]])) {
    export(smoothed_ls_bigwig[[i]][[j]], 
           paste0(output_path, 
                  smoothing_window_ls[i], 
                  "bp/bigwig/Rec12_IP_vs_INPUT_", 
                  names(smoothed_ls_bigwig[[i]])[j], 
                  "_", 
                  smoothing_window_ls[i], 
                  "bp_smoothed.bw"), 
           format = "BigWig")
  }
}



