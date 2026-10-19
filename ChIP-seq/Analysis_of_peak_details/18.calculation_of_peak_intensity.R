#function to calculate the intensity (and density) of ChIP-seq in the defined regions
#results are indicated in "Table_S5_NGS_data_in_hotspots.xlsx"

#datasets before smoothing are used 
#before calculation, appropriate mask data (defined in 11.mask_based_on_INPUT.R) are applied.

library(dplyr)
library(data.table)

peak_intensity <- function(peak_data, coverage_data_ls, cal_density = FALSE){
  
  #cal peak intensity
  peak_intensity <- lapply(1:nrow(peak_data), function(n){
    current_chr         <- peak_data$chr[n]
    current_start       <- peak_data$start[n]
    current_end         <- peak_data$end[n]
    current_peak_length <- current_end - current_start
    
    count_in_peak <- coverage_data_ls[[current_chr]]$count[current_start:current_end] %>% 
      as.numeric()
    
    total_count_in_peak <- count_in_peak %>% 
      sum(na.rm = TRUE)
    
    if(isTRUE(cal_density)){
      density <- total_count_in_peak / current_peak_length
      result_df <- data.frame(
        name    = peak_data$name[n],
        density = density
      ) 
    }else{
      result_df <- data.frame(
        name      = peak_data$name[n],
        intensity = total_count_in_peak
      ) 
    }
    
    return(result_df)  
  })
  
  peak_intensity_df <- do.call(rbind, peak_intensity)
  return(peak_intensity_df)
}


#e.g.,
#hotspot data
peak_df <- fread("WT_glu_plus_minus_rec12_merged.bed")

#coverage data
load("Rec12_ChIP_seq_IP_vs_INPUT_10bp_bin_masked.RData.gz")
coverage_ls <- IP_vs_INPUT_normalized_peak_SK1_ls

#calculate
Rec12_intensity_hotspot <- lapply(coverage_ls, function(n){
  peak_intensity(peak_data        = peak_df,
                 coverage_data_ls = n,
                 cal_density      = FALSE)
})