setwd("")
library(data.table)
library(dplyr)
library(rtracklayer)

#data
base_path <- ""
path <- paste0(base_path, "/coverage_data/masked_10bp_binned/Rdata/")
output_path <- paste0(base_path, "/coverage_data/normalized/")

output_name <- c(
  ""
)

INPUT_index <- 1:4
IP_index    <- 5:8


load(paste0(path, "/Rec12_ChIP_seq_coverage_ls_10bp_bin_masked.RData.gz"))

#normalization by spike-in
spike_in_SK1_in_peaks <- sapply(names(masked_data_binned_10bp), function(n){
  SK1_read <- fread(paste0(base_path,
                           "/bam_file/bam_SK1/SK1_peaks/SK1_peak_read/mapping_status/",
                           n, "_mapped_status_SK1_peaks.txt"))
  SK1_read_wo_chrMT_peaks <- SK1_read %>% 
    filter(grepl("SK1_", V1)) %>% 
    filter(V1 != "SK1_MT") %>% 
    pull(V3) %>% 
    sum()
})


#normalize
normalized_data_ls <- lapply(1:length(masked_data_binned_10bp), function(n){
  result <- lapply(masked_data_binned_10bp[[n]], function(i){
    normalized_data <- i$count / spike_in_SK1_in_peaks[n] * 1000000
    result_df       <- i
    result_df$count <- normalized_data
    return(result_df)
  })
  names(result) <- names(masked_data_binned_10bp[[n]])
  return(result)
})
names(normalized_data_ls) <- names(masked_data_binned_10bp)

#cal IP/INPUT
#N/0 is replaced by NA

IP_vs_INPUT_data_ls <- lapply(1:length(IP_index), function(n){
  INPUT_data <- normalized_data_ls[[INPUT_index[n]]]
  IP_data    <- normalized_data_ls[[IP_index[n]]]
  
  result <- lapply(1:length(INPUT_data), function(i){
    IP_vs_INPUT <- IP_data[[i]]$count / INPUT_data[[i]]$count
    IP_vs_INPUT[is.nan(IP_vs_INPUT)]      <- NA
    IP_vs_INPUT[is.infinite(IP_vs_INPUT)] <- NA
    result_df       <- IP_data[[i]]
    result_df$count <- IP_vs_INPUT
    return(result_df)
  })
  names(result) <- names(INPUT_data)
  return(result)
})
names(IP_vs_INPUT_data_ls) <- output_name


#save as RData
normalized_peak_SK1_ls <- normalized_data_ls
IP_vs_INPUT_normalized_peak_SK1_ls <- IP_vs_INPUT_data_ls

dir.create(output_path)
dir.create(paste0(output_path, "/RData/"))
dir.create(paste0(output_path, "/RData/coverage/"))
dir.create(paste0(output_path, "/RData/IP_vs_INPUT/"))

save(normalized_peak_SK1_ls,             file = paste0(output_path, "/RData/coverage/Rec12_ChIP_seq_coverage_ls_10bp_bin_masked_normalized_peak_SK1.RData.gz"), compress = TRUE)
save(IP_vs_INPUT_normalized_peak_SK1_ls, file = paste0(output_path, "/RData/IP_vs_INPUT/Rec12_ChIP_seq_IP_vs_INPUT_10bp_bin_masked.RData.gz"), compress = TRUE)

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
normalized_peak_SK1_ls_bigwig             <- lapply(normalized_peak_SK1_ls, function(n){
  change_format(data_ls  = n, 
                used_col = "count"
  )
})
IP_vs_INPUT_normalized_peak_SK1_ls_bigwig <- lapply(IP_vs_INPUT_normalized_peak_SK1_ls, function(n){
  change_format(data_ls  = n, 
                used_col = "count"
  )
})

#save as bigwig
dir.create(paste0(output_path, "/bigwig/"))
dir.create(paste0(output_path, "/bigwig/coverage/"))
dir.create(paste0(output_path, "/bigwig/IP_vs_INPUT/"))

for (i in 1:length(normalized_peak_SK1_ls_bigwig)) {
  export(normalized_peak_SK1_ls_bigwig[[i]], 
         paste0(output_path, "/bigwig/coverage/", names(normalized_peak_SK1_ls_bigwig)[i], ".bw"), 
         format = "BigWig")
}

for (i in 1:length(IP_vs_INPUT_normalized_peak_SK1_ls_bigwig)) {
  export(IP_vs_INPUT_normalized_peak_SK1_ls_bigwig[[i]], 
         paste0(output_path, "/bigwig/IP_vs_INPUT/", names(IP_vs_INPUT_normalized_peak_SK1_ls_bigwig)[i], ".bw"), 
         format = "BigWig")
}
