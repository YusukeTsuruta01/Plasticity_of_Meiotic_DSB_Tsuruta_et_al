setwd("")
library(dplyr)
library(data.table)

window_size <- 501

base_path <- ""
path <- paste0(base_path, "/coverage_data/average/", window_size, "bp/RData/")
output_path <- paste0(base_path, "/coverage_data/average/")

output_name <- c(
  ""
)

#average data
#501bp smoothed
load(paste0(path, 
            "/Rec12_IP_vs_INPUT_",
            window_size, 
            "_smoothed_average.RData.gz"))


#average of IP/INPUT
ave <- sapply(Rec12_IP_vs_INPUT_smoothed_ave, function(n){
  mean(do.call(rbind, n)$count, na.rm = TRUE)
})
#use average x1.5 as a threshold

#function for peak calling
peak_calling <- function(data, peak_threshold, chr_col = "chr", score_col = "count"){
  #data
  score           <- pull(data, !!sym(score_col))
  current_chr     <- pull(data, !!sym(chr_col))[1]
  
  #mark index which is more than threshold
  over_threshold <- ifelse(score >= peak_threshold, 1, 0)
  over_threshold[is.na(over_threshold)] <- 0
  
  #define peak edge
  over_threshold_lag  <- lag(over_threshold,  n = 1, default = 0)
  peak_edge_start <- which((over_threshold - over_threshold_lag) == 1)
  
  over_threshold_lead <- lead(over_threshold, n = 1, default = 0)
  peak_edge_end  <- which((over_threshold - over_threshold_lead) == 1)
  
  #define peak summit and trough
  #code is from Dr.Murakami
  peaksign <- function(series,span=3,do.pad=TRUE) {
    # Purpose: return (-1 / 0 / 1) if series[i] is ( trough / "normal" / peak )
    # Author: Martin Maechler, Date: 25 Nov 2005
    if((span<-as.integer(span))%%2!=1 || span==1) stop("'span' must be odd and >= 3")
    s1=1+(s<-span%/%2)
    z=embed(series,span)
    d=z[,s1]-z[,-s1]
    ans=rep.int(0,nrow(d))
    ans[apply(d>=0,1,all)]=as.integer(1)
    ans[apply(d<=0,1,all)]=as.integer(-1)
    if(do.pad) {
      pad=rep.int(0,s)
      c(pad,ans,pad)
    } else ans
  }
  ##
  
  #peak sumittes are 1, troughs are -1, egdes are -2, others are 0
  peak_summit_and_trough <- peaksign(series = score, span = 3, do.pad = TRUE)
  peak_summit_and_trough[which(over_threshold == 0)] <- 0
  peak_summit_and_trough[peak_edge_start] <- -2
  peak_summit_and_trough[peak_edge_end]   <- -2
  
  #define peak region
  index_trough_and_egde <- which(peak_summit_and_trough %in% c(-1, -2))
  if(length(index_trough_and_egde) > 0){
    result_ls     <- lapply(1:(length(index_trough_and_egde)-1), function(n){
      
      current_index <- index_trough_and_egde[n]
      next_index    <- index_trough_and_egde[n + 1]
      
      current_peak_region <- peak_summit_and_trough[current_index:next_index]
      
      if(any(current_peak_region == 1)){
        result_df <- data.frame(
          start   = current_index,
          end     = next_index,
          sumitte = current_index + which(current_peak_region == 1) - 1
        )
      }else{
        result_df <- data.frame(
          start   = integer(0),
          end     = integer(0),
          sumitte = integer(0)
        )
      }
      return(result_df)
    })
    
    #data.frame
    final_peak_df   <- do.call(rbind, result_ls)
    final_peak_df_2 <- final_peak_df %>% 
      mutate(chr = rep(current_chr, times = nrow(final_peak_df))) %>% 
      select(chr, start, end, sumitte)
  }else{
    final_peak_df_2 <- data.frame(
      start   = integer(0),
      end     = integer(0),
      sumitte = integer(0)
    )
  }
  return(final_peak_df_2)
}

peak_ls <- lapply(1:length(Rec12_IP_vs_INPUT_smoothed_ave), function(n){
  
  result_ls <- lapply(1:length(Rec12_IP_vs_INPUT_smoothed_ave[[n]]), function(i){
    peak_calling(data = Rec12_IP_vs_INPUT_smoothed_ave[[n]][[i]], 
                 peak_threshold = ave[n] * 1.5, 
                 chr_col = "chr", 
                 score_col = "count")
  })
  
  result_ls_df   <- do.call(rbind, result_ls)
  result_ls_df_2 <- result_ls_df %>% 
    mutate(name = paste0(
      names(Rec12_IP_vs_INPUT_smoothed_ave)[n], "_rec12_peaks_", 1:nrow(result_ls_df)
    ))
  
  return(result_ls_df_2)
})
names(peak_ls) <- output_name

#make bed file to check by IGV
dir.create(paste0(base_path, "/analysis/"))
dir.create(paste0(base_path, "/analysis/peak_call/"))
dir.create(paste0(base_path, "/analysis/peak_call/before_merge/"))

for (i in 1:length(peak_ls)) {
  output_file <- NULL
  output_file <- peak_ls[[i]] %>% 
    select(chr, start, end, name)
  write.table(output_file, 
              file = paste0(base_path, 
                            "/analysis/peak_call/before_merge/",
                            names(peak_ls)[i], 
                            "_before_merge.bed"), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}


