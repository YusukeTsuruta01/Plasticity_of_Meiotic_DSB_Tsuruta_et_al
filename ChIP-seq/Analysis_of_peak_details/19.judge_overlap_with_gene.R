#function to judge whether hotspots are overlapped with genes
#results are indicated in "Table_S5_NGS_data_in_hotspots.xlsx"

library(dplyr)
library(data.table)

judge_overlap_with_RNA <- function(peak_data, RNA_data, 
                                   distance_threshold = 1500, 
                                   gene_type = "protein_coding_gene", 
                                   name_col  = 4,
                                   RNA_seq_data_col = 11:12){
  #data of RNA-seq
  filtered_RNA_data <- RNA_data %>% 
    filter(type == gene_type)
  
  #coordinate of upstream
  upstream_start <- ifelse(filtered_RNA_data$strand == "+", 
                           filtered_RNA_data$start - distance_threshold,
                           filtered_RNA_data$end   + 1)
  
  upstream_end <- ifelse(filtered_RNA_data$strand == "+", 
                         filtered_RNA_data$start - 1,
                         filtered_RNA_data$end   + distance_threshold)
  
  filtered_RNA_data_upstream <- filtered_RNA_data %>% 
    mutate(gene_id = paste0(gene_id, "_upstream"),
           start   = upstream_start,
           end     = upstream_end,
    )
  
  filtered_RNA_data_2 <- rbind(filtered_RNA_data, filtered_RNA_data_upstream)
  
  #data of peaks
  peak_chr   <- peak_data$chr
  peak_name  <- peak_data[[name_col]]
  peak_start <- peak_data$start
  peak_end   <- peak_data$end
  
  
  result <- lapply(1:nrow(peak_data), function(n){
    crrent_peak_chr   <- peak_chr[n]
    crrent_peak_start <- peak_start[n]
    crrent_peak_end   <- peak_end[n]
    
    candidate <- filtered_RNA_data_2 %>% 
      filter(chr == crrent_peak_chr) %>% 
      filter(end >= crrent_peak_start) %>% 
      filter(start <= crrent_peak_end)
    
    if(nrow(candidate) == 0){
      candidate <- candidate %>% 
        select(1, all_of(RNA_seq_data_col))
      #print(candidate)
      empty <- data.frame(matrix(NA, nrow = 1, ncol = ncol(candidate)))
      colnames(empty) <- colnames(candidate)
      candidate <- empty
    }else if(nrow(candidate) >= 1){
      overlap_length_ls <- sapply(1:nrow(candidate), function(i){
        coordinate_ls <- c(crrent_peak_start, crrent_peak_end, candidate$end[i], candidate$start[i]) %>% 
          sort()
        overlap_length <- coordinate_ls[3] - coordinate_ls[2]
      })
      overlap_rate <- overlap_length_ls / (crrent_peak_end - crrent_peak_start + 1)
      
      if(max(overlap_rate) < 0.25){
        candidate <- candidate %>% 
          select(1, all_of(RNA_seq_data_col))
        #print(candidate)
        empty <- data.frame(matrix(NA, nrow = 1, ncol = ncol(candidate)))
        colnames(empty) <- colnames(candidate)
        candidate <- empty
      }else{
        #print(candidate$gene_id)
        #print(overlap_rate)
        max_index <- which(overlap_rate == max(overlap_rate))
        if(length(max_index) == 1){
          candidate <- candidate[max_index, ] %>% 
            select(1, all_of(RNA_seq_data_col))
        }else{
          candidate <- candidate[max_index, ] %>% 
            filter(TPM_glu_minus == max(TPM_glu_minus, na.rm = TRUE)) %>% 
            select(1, all_of(RNA_seq_data_col))
        }
        #print(candidate) 
      }
    }
    
    candidate <- candidate %>% 
      mutate(name = peak_name[n])
    #print(candidate)
    return(candidate)
  }) 
  #print(result)
  #print(sapply(result, function(n){ncol(n)}))
  result_df <- do.call(rbind, result)
  
  result_df <- result_df %>% 
    select(name, everything())
  
  return(result_df)
}

#e.g.,
#judge overlap with hotspot and protein coding gene
peak_df <- fread("WT_glu_plus_minus_rec12_merged.bed")
RNA_seq <- fread("RNA_seq_WT_meiosis_glu_plus_minus_N1_260608.txt")
#This is the result of feature counts normalized by TPM
#The coordinate in the result of feature counts contain UTR

closest_protein_coding_gene <- judge_overlap_with_RNA(peak_data = peak_df, 
                                                      RNA_data  = pombe_RNAseq, 
                                                      distance_threshold = 1500, 
                                                      gene_type = "protein_coding_gene",
                                                      name_col  = 4,
                                                      RNA_seq_data_col = 11:12)
