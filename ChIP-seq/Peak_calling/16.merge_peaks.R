setwd("")
library(dplyr)
library(data.table)

base_path <- ""
output_dr <- paste0(base_path, "/analysis/peak_call/merged_peak/")


peak_file_name <- c(
  ""
)

peak_df_plus  <- fread(paste0(base_path, "/analysis/peak_call/before_merge/", peak_file_name[1], "_before_merge.bed"))
peak_df_minus <- fread(paste0(base_path, "/analysis/peak_call/before_merge/", peak_file_name[2], "_before_merge.bed"))

#merge peaks between plus and minus
#retaining sub-peaks

#check overlap
check_overlap <- function(peak_df_1, peak_df_2, overlap_ratio){
  final_result <- lapply(1:nrow(peak_df_1), function(n){
    
    current_chr   <- peak_df_1$V1[n]
    current_start <- peak_df_1$V2[n]
    current_end   <- peak_df_1$V3[n]
    current_name  <- peak_df_1$V3[n]
    
    candidate <- peak_df_2 %>% 
      filter(V1 == current_chr) %>% 
      filter(V2 <= current_end & V3 >= current_start)
    
    if(nrow(candidate) != 0){
      coordinate_df  <- candidate %>% 
        mutate(current_start = rep(current_start, times = nrow(candidate)),
               current_end   = rep(current_end  , times = nrow(candidate)),
               current_name  = rep(current_name  , times = nrow(candidate)))
      
      overlapped_peak <- lapply(1:nrow(candidate), function(i){
        coordinate <- c(coordinate_df$V2[i], 
                        coordinate_df$V3[i],
                        coordinate_df$current_start[i],
                        coordinate_df$current_end[i]) %>% sort()
        
        overlap_length <- coordinate[3] - coordinate[2]
        
        overlap_ratio_df1 <- overlap_length / (coordinate_df$current_end[i] - coordinate_df$current_start[i])
        overlap_ratio_df2 <- overlap_length / (coordinate_df$V3[i] - coordinate_df$V2[i])
        
        if(overlap_ratio_df1 > overlap_ratio | overlap_ratio_df2 > overlap_ratio){
          result <- peak_df_1[n, ] %>% 
            mutate(
              overlapped_by         = candidate$V4[i],
              overlapped_peak_start = candidate$V2[i],
              overlapped_peak_end   = candidate$V3[i]
            )
          names(result) <- c("chr", "start", "end", "names", "overlapped_by", "overlapped_peak_start", "overlapped_peak_end")
          #print(result)
        }else{
          result <-  peak_df_1[n, ] %>% 
            mutate(
              overlapped_by         = as.character(NA),
              overlapped_peak_start = NA,
              overlapped_peak_end   = NA
            )
          names(result) <- c("chr", "start", "end", "names", "overlapped_by", "overlapped_peak_start", "overlapped_peak_end")
        }
        return(result)
      })
      overlapped_peak_df <- do.call(rbind, overlapped_peak)
      #print(overlapped_peak_df)
    }else{
      overlapped_peak_df <- peak_df_1[n, ] %>% 
        mutate(
          overlapped_by         = as.character(NA),
          overlapped_peak_start = NA,
          overlapped_peak_end   = NA
        )
      names(overlapped_peak_df) <- c("chr", "start", "end", "names", "overlapped_by", "overlapped_peak_start", "overlapped_peak_end")
    }
    return(overlapped_peak_df)
    #print(overlapped_peak_df)
  })
  final_result_df <- do.call(rbind, final_result)
  #print(final_result_df %>% filter(names == "WT_glu_plus_rec12_peaks_396"))
  
  final_result_wo_dup <- lapply(unique(final_result_df$names), function(n){
    result_with_dup <- final_result_df %>% 
      filter(names == n)
    
    if(nrow(result_with_dup) == 1){
      result <- result_with_dup
    }else if(all(is.na(result_with_dup$overlapped_by))){
      result <- result_with_dup[1, ]
    }else{
      result <- result_with_dup %>% 
        filter(!is.na(overlapped_by))
    }
    return(result)
  })
  final_result_wo_dup_df <- do.call(rbind, final_result_wo_dup)
  
  return(final_result_wo_dup_df)
}

plus_overlapped_by_minus <- check_overlap(peak_df_plus, 
                                          peak_df_minus, 
                                          overlap_ratio = 0.5)
minus_overlapped_by_plus <- check_overlap(peak_df_minus, 
                                          peak_df_plus, 
                                          overlap_ratio = 0.5)

merge_overlapped_peaks <- function(peak_df1, peak_df2){
  #no overlap
  no_overlap_peak_df1 <- peak_df1 %>% 
    filter(is.na(overlapped_by)) %>% 
    select(chr, start, end, names)
  no_overlap_peak_df2 <- peak_df2 %>% 
    filter(is.na(overlapped_by)) %>% 
    select(chr, start, end, names)
  
  #overlap(1:1)
  overlap_1by1 <- lapply(unique(peak_df1$names), function(n){
    candidata <- peak_df1 %>% 
      filter(names == n) %>% 
      filter(!is.na(overlapped_by))
    #print(candidata)
    
    if(nrow(candidata) == 1){
      if((nrow(peak_df2[peak_df2$names == candidata$overlapped_by])) == 1){
        coordinate <- c(candidata$start, candidata$end,
                        candidata$overlapped_peak_start,
                        candidata$overlapped_peak_end)
        #print(coordinate)
        overlap_1by1 <- data.frame(
          chr   = candidata$chr,
          start = min(coordinate),
          end   = max(coordinate),
          names = paste0(candidata$names, ",", candidata$overlapped_by)
        )
      }else{
        overlap_1by1 <- data.frame(
          chr   = as.character(NA),
          start = NA,
          end   = NA,
          names = as.character(NA)
        )
      }
    }else{
      overlap_1by1 <- data.frame(
        chr   = as.character(NA),
        start = NA,
        end   = NA,
        names = as.character(NA)
      )
    }
    return(overlap_1by1)
  })
  overlap_1by1_df <- do.call(rbind, overlap_1by1) %>% 
    filter(!is.na(names))
  
  #overlap(1:2~)
  total_df <- rbind(peak_df1, peak_df2)
  overlap_1by_more_than_2_candidate <- lapply(unique(total_df$names), function(n){
    candidate_df <- total_df %>% 
      select(chr, start, end, names, overlapped_by) %>% 
      filter(!is.na(overlapped_by)) %>% 
      filter(names == n)
    
    if(nrow(candidate_df) > 1){
      return(candidate_df)
    }else{
      candidate_df <- data.frame(
        chr   = as.character(NA),
        start = NA,
        end   = NA,
        names = as.character(NA),
        overlapped_by = as.character(NA)
      )
      return(candidate_df)
    }
  })
  overlap_1by_more_than_2_candidate_df <- do.call(rbind, overlap_1by_more_than_2_candidate) %>% 
    filter(!is.na(names))
  
  #grouping of overlapping peaks
  overlapping_peak_group <- lapply(unique(overlap_1by_more_than_2_candidate_df$names), function(n){
    current_peak <- n
    result_vec   <- n
    
    repeat {
      next_peak <- overlap_1by_more_than_2_candidate_df %>% 
        filter(names == current_peak) %>% 
        pull(overlapped_by)
      
      if (length(next_peak) == 0 || all(next_peak %in% result_vec)) break
      
      result_vec   <- c(result_vec, next_peak)
      result_vec   <- unique(result_vec)
      current_peak <- next_peak[length(next_peak)]
    }
    return(result_vec)
  })
  
  #remove duplicate
  overlapping_peak_group_wo_dup <- lapply(1:length(overlapping_peak_group), function(n){
    current_group <- overlapping_peak_group[[n]]
    terget_group  <- overlapping_peak_group[-n]
    
    chech_dup <- sapply(terget_group, function(i){
      all(current_group %in% i)
    })
    if(!any(chech_dup)){
      result <- current_group
    }else{
      result <- NA
    }
  })
  overlapping_peak_group_wo_dup <- overlapping_peak_group_wo_dup[!is.na(overlapping_peak_group_wo_dup)]
  
  #make new sub peak region
  sub_peak_define <- lapply(overlapping_peak_group_wo_dup, function(n){
    #plus peak
    plus_peaks         <- grep("plus",  n, value = TRUE)
    plus_peaks_detail  <- peak_df1 %>% 
      filter(names %in% plus_peaks) %>% 
      select(chr, start, end, names) %>% 
      unique()
    #print(plus_peaks_detail)
    
    if(length(plus_peaks) >= 2){
      plus_peaks_junction <- sapply(1:(length(plus_peaks) -1), function(n){
        junction <- rep(
          round((plus_peaks_detail$end[n] + plus_peaks_detail$start[n+1])/2),
          times = 2)
      }) 
    }else{
      plus_peaks_junction <- NA
    }
    
    
    #minus peak
    minus_peaks         <- grep("minus",  n, value = TRUE)
    minus_peaks_detail  <- peak_df2 %>% 
      filter(names %in% minus_peaks) %>% 
      select(chr, start, end, names) %>% 
      unique()
    #print(minus_peaks_detail)
    if(length(minus_peaks) >= 2){
      minus_peaks_junction <- sapply(1:(length(minus_peaks) -1), function(n){
        junction <- rep(
          round((minus_peaks_detail$end[n] + minus_peaks_detail$start[n+1])/2),
          times = 2)
      }) 
    }else{
      minus_peaks_junction <- NA
    }
    
    #print(plus_peaks_detail)
    #print(minus_peaks_detail)
    #print(min(c(plus_peaks_detail$start, minus_peaks_detail$start)))
    #print(plus_peaks_junction)
    
    #total coordinate with junction
    total_coordinate <- c(
      min(c(plus_peaks_detail$start, minus_peaks_detail$start)),
      max(c(plus_peaks_detail$end,   minus_peaks_detail$end)),
      plus_peaks_junction,
      minus_peaks_junction
    )
    
    total_coordinate <- total_coordinate[!is.na(total_coordinate)] %>% 
      sort()
    #print(total_coordinate)
    total_name <- paste0(n, collapse = ",")
    
    coordinate_df <- matrix(total_coordinate, 
                            ncol = 2, 
                            byrow = TRUE) %>% 
      as.data.frame()
    names(coordinate_df) <- c("start", "end")
    
    result <- coordinate_df %>% 
      mutate(
        chr   = rep(unique(plus_peaks_detail$chr), times = nrow(coordinate_df)),
        names = rep(total_name,                    times = nrow(coordinate_df))
      ) %>% 
      select(chr, start, end, names)
    return(result)
  })
  sub_peak_define_df <- do.call(rbind, sub_peak_define)
  
  #use only chrI, II, III
  final_result <- rbind(no_overlap_peak_df1,
                        no_overlap_peak_df2,
                        overlap_1by1_df,
                        sub_peak_define_df) %>% 
    filter(chr %in% c("I", "II", "III"))
  
  #sort
  final_result$chr    <- factor(final_result$chr, level = c("I", "II", "III"))
  final_result_sorted <- final_result %>% 
    arrange(chr, start) %>% 
    mutate(names2 = paste0("rec12_peaks_", 1:nrow(final_result))) %>% 
    select(chr, start, end, names2, names)
  final_result_sorted$chr <- as.character(final_result_sorted$chr)
  names(final_result_sorted) <- c("chr", "start", "end", "names", "overlapped")
  
  return(final_result_sorted)
}

merged_peak <- merge_overlapped_peaks(plus_overlapped_by_minus, minus_overlapped_by_plus)

dir.create(output_dr)
write.table(merged_peak, file = paste0(output_dr, "/WT_glu_plus_minus_rec12_merged.bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
