setwd("")
library(data.table)
library(dplyr)
library(rtracklayer)

#data
base_path <- ""
path <- paste0(base_path, "/coverage_data/10bp_binned/Rdata/")
output_path <- paste0(base_path, "/coverage_data/masked_10bp_binned/")

load(paste0(path, "/Rec12_ChIP_seq_coverage_ls_10bp_bin.RData.gz"))

INPUT_index <- 1:4
chromosome  <- c("I", "II", "III")


INPUT_mean  <- sapply(binned_10bp[INPUT_index], function(n){
  current_data <- do.call(rbind, n) %>% 
    filter(chr %in% chromosome) %>% 
    pull(count)
  mean(current_data)
})

INPUT_sd  <- sapply(binned_10bp[INPUT_index], function(n){
  current_data <- do.call(rbind, n) %>% 
    filter(chr %in% chromosome) %>% 
    pull(count)
  sd(current_data)
})
#INPUT_mean +/- 2.5*INPUT_sd might be good for WT rep1,2

#function for mask (from Dr. Murakami)
get.mask = function(dat, mincov, maxcov, extender=1000){
  mydat <- double(length(dat))
  mydat[1:length(dat)] <- dat
  tmp <- ifelse(mydat<mincov | mydat>maxcov, 0, 1)
  #leftbounds <- c(1, which(diff(tmp)==-1))
  if(tmp[1]==1) {leftbounds = c(which(diff(tmp)==-1))} else {leftbounds = c(1, which(diff(tmp)==-1))}
  
  #rightbounds <- c(which(diff(tmp)==1), length(tmp))
  if(tmp[length(tmp)]==1) {rightbounds = c(which(diff(tmp)==1))} else {rightbounds = c(which(diff(tmp)==1), length(tmp))}
  leftbounds <- ifelse(leftbounds<extender,1,leftbounds-extender)
  rightbounds <- ifelse(rightbounds>(length(tmp)-extender), length(tmp), rightbounds+extender)
  for(i in 1:length(leftbounds)){
    tmp[leftbounds[i]:rightbounds[i]] <- 0
  }
  inds = which((rightbounds-leftbounds+1) <= (2*extender+10)) ### this line was originally design to filter out masks less than 2*extender + 10 bp but this has problem. anyway modification of assiging 0 inactivated this. For the future, if the activation is needed, this needs to be done before extending masks.
  if(length(inds)>0){
    for(i in inds) tmp[leftbounds[i]:rightbounds[i]] <- 0
  }
  tmp <- ifelse(tmp==1,1,NA)
  return(tmp)
}

#mask
mask_region_ls <- lapply(INPUT_index, function(n){
  current_data <- do.call(rbind, binned_10bp[[n]])
  get.mask(dat = current_data$count, 
           mincov = INPUT_mean[n] - INPUT_sd[n] * 2.5, 
           maxcov = INPUT_mean[n] + INPUT_sd[n] * 2.5, 
           extender = 250
           )
})
merged_mask_region <- Reduce("*", mask_region_ls)
chr_index <- do.call(rbind, binned_10bp[[1]])$chr
merged_mask_region_ls <- lapply(unique(chr_index), function(n){
  merged_mask_region[chr_index == n]
})
merged_mask_region_ls[[4]] <- rep(NA, times = 20000)
merged_mask_region_ls[[5]] <- rep(NA, times = 19431)
merged_mask_region_ls[[6]] <- rep(NA, times = 20128)

#save mask
dir.create(output_path)
save(merged_mask_region_ls, file = paste0(output_path, "/mask.RData"))

#mask data
masked_data_binned_10bp <- lapply(binned_10bp, function(n){
  result_ls <- lapply(1:length(n), function(i){
    result <- n[[i]]
    result$count <- result$count * merged_mask_region_ls[[i]]
    return(result)
  })
  names(result_ls) <- names(n)
  return(result_ls)
})

#save as RData
dir.create(paste0(output_path, "/RData/"))
save(masked_data_binned_10bp, file = paste0(output_path, "/RData/Rec12_ChIP_seq_coverage_ls_10bp_bin_masked.RData.gz"), compress = TRUE)

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
masked_data_binned_10bp_bigwig <- lapply(masked_data_binned_10bp, function(n){
  change_format(data_ls  = n, 
                used_col = "count"
  )
})

#save as bigwig
dir.create(paste0(output_path, "/bigwig/"))
for (i in 1:length(masked_data_binned_10bp_bigwig)) {
  export(masked_data_binned_10bp_bigwig[[i]], 
         paste0(output_path, "/bigwig/", names(masked_data_binned_10bp_bigwig)[i], ".bw"), 
         format = "BigWig")
}

