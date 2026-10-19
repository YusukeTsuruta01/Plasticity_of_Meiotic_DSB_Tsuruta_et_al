setwd("")
library(data.table)
library(dplyr)
library(rtracklayer)

#data
base_path <- ""
path <- paste0(base_path, "/coverage_data/from_bam_file/")
output_path <- paste0(base_path, "/coverage_data/10bp_binned/")

sample_name <- c(
  ""
)

raw_data_ls <- lapply(sample_name, function(n){
  data <- fread(paste0(path, n, "_coverage_map.txt.gz"))
})
names(raw_data_ls) <- sample_name
chr_name <- unique(raw_data_ls[[1]]$V1)[1:6]

#binning before calculate IP vs INPUT
binnig <- function(df, window, chr_length){
  
  bin_number <- ifelse(chr_length %% window == 0, chr_length %/% window, chr_length %/% window + 1)
  
  bin_start  <- 0:(bin_number-1) *window + 1
  bin_end    <- 1:bin_number *window
  bin_end[length(bin_end)] <- chr_length
  
  bin_ave <- sapply(1:bin_number, function(n){
    mean(df$V3[bin_start[n]:bin_end[n]], na.rm = TRUE)
  })
  bin_ave[is.nan(bin_ave)] <- NA
  
  last_bin_size <- bin_end[length(bin_end)] - bin_start[length(bin_start)] + 1
  
  result_df <- data.frame(
    chr        = rep(df$V1[1], times = chr_length),
    coordinate = 1:chr_length,
    count      = c(
      rep(bin_ave[1:(bin_number - 1)], each  = window), 
      rep(bin_ave[bin_number],         times = last_bin_size)
    )
  )
  return(result_df)
}
binned_10bp <- lapply(raw_data_ls, function(n){
  result <- lapply(chr_name, function(i){
    current_data       <- n[n$V1 == i, ]
    current_chr_length <- nrow(current_data)
    
    result_df <- binnig(df = current_data, 
                        window = 10, 
                        chr_length = current_chr_length)
    return(result_df)
  })
  names(result) <- chr_name
  return(result)
})

#save as RData
dir.create(output_path)
dir.create(paste0(output_path, "/RData/"))
save(binned_10bp, file = paste0(output_path, "/RData/Rec12_ChIP_seq_coverage_ls_10bp_bin.RData.gz"), compress = TRUE)

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
binned_10bp_bigwig <- lapply(binned_10bp, function(n){
  change_format(data_ls  = n, 
                used_col = "count"
                )
})

#save as bigwig
dir.create(paste0(output_path, "/bigwig/"))
for (i in 1:length(binned_10bp_bigwig)) {
  export(binned_10bp_bigwig[[i]], 
         paste0(output_path, "/bigwig/", names(binned_10bp_bigwig)[i], ".bw"), 
         format = "BigWig")
}



