setwd("")
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(rtracklayer)
library(zoo)

#Line plot to show the ChIP-seq profile

#This code is used for the following figure
#Fig. 1C, 1D, S1B, S1C,, S1F
#Fig. 3A
#Fig. 4A, S4E 
#Fig. 5B, 5E, S5D, S5E, S5F
#Fig. 6C, S6E, S6F, S6G, S6H, S6I

#function to show the ChIP-seq profile
profile_plot <- function(count_data,
                         count_data_col,
                         chromosome,
                         sample_name_ls,
                         xlim   = NULL, 
                         ylim   = NULL, 
                         x_label = NULL, 
                         y_label = NULL, 
                         title_name = NULL, 
                         binned = FALSE,
                         na_approx = FALSE,
                         font_size   = 12,
                         plot_margin = c(0,0,0,0)){
  
  if(isTRUE(binned)){
    count_data <- lapply(count_data, function(n){
      n %>% 
        mutate(coordinate = (bin_start + bin_end) / 2)
    })
  }
  
  x_coordinate <- count_data[[1]] %>% 
    filter(chr == chromosome) %>% 
    filter(coordinate >= min(xlim) & coordinate <= max(xlim)) %>% 
    pull(coordinate)
  
  #count data
  count_data_ls <- lapply(count_data, function(n){
    plot_data <- n %>% 
      filter(chr == chromosome) %>% 
      filter(coordinate >= min(xlim) & coordinate <= max(xlim)) %>% 
      pull(!!sym(count_data_col))
  })
  
  if(isTRUE(na_approx)){
    count_data_ls <- lapply(count_data, function(n){
      plot_data <- n %>% 
        filter(chr == chromosome)
      
      plot_data[[count_data_col]][plot_data[[count_data_col]] == 0] <- NA
      plot_data[[count_data_col]] <- na.approx(plot_data[[count_data_col]], na.rm = FALSE)
      
      plot_data <- plot_data %>% 
        filter(coordinate >= min(xlim) & coordinate <= max(xlim)) %>% 
        pull(!!sym(count_data_col))
    })
  }
  
  #df for plot
  plot_data <- do.call(cbind, count_data_ls) %>% 
    as.data.table() %>% 
    mutate(x = x_coordinate)
  names(plot_data) <- c(sample_name_ls, "x")
  
  plot_data <- data.table::melt(plot_data, 
                                id.vars = "x",
                                variable.name = "sample", 
                                value.name = "count")
  
  plot_data$sample <- factor(plot_data$sample, levels = sample_name_ls)
  
  final_fig <- ggplot()+
    geom_line(data = plot_data, aes(x = x, y = count, colour = sample), linewidth = 0.75)+
    scale_x_continuous(limits = xlim, expand = c(0,0), labels = scales::label_number(accuracy = 1))+
    scale_y_continuous(limits = ylim, expand = c(0,0), labels = scales::label_number(accuracy = 0.1))+
    labs(x = x_label, 
         y = y_label, 
         colour = NULL, 
         title = title_name)+
    theme_classic(base_size = font_size)+
    theme(
      text        = element_text(size = font_size, color = "black"),
      axis.text   = element_blank(),
      legend.text = element_text(size = font_size, color = "black"),
      plot.title  = element_text(size = font_size, color = "black", hjust = 0.5),
      legend.position = "none",
      plot.margin = unit(plot_margin, "cm")
    )
  return(final_fig)
}

#function to show the gene
gene_plot <- function(gene, 
                      chr, 
                      xlim, 
                      font_size = 12,
                      plot_margin = c(0,0,0,0)){
  
  #V1, chromosome
  #V2, start
  #V3, end
  #V5, strand
  
  current_data <- gene %>% 
    filter(V1 == chr)
  
  ggplot()+
    geom_segment(aes(x = min(xlim), xend = max(xlim), 
                     y = 0, yend = 0))+
    geom_rect(data = current_data, 
              aes(xmin = V2, xmax = V3, 
                  ymin = -1, ymax = 1, fill = V6))+
    geom_text(data = current_data, 
              aes(x = V2,label = V4, angle = 90), y = 2)+
    scale_x_continuous(limits = c(min(xlim), max(xlim)), expand = c(0,0))+
    scale_y_continuous(limits = c(-1, 3))+
    theme_classic()+
    theme(
      text        = element_text(size = font_size, color = "black"),
      axis.text   = element_blank(),
      legend.text = element_text(size = font_size, color = "black"),
      plot.title  = element_blank(),
      legend.position = "none",
      plot.margin = unit(plot_margin, "cm")
    )
}


#e.g.,
#Fig. 1B

##Processing data
sample_name <- c(
  "GSM10129174_Rec12_WT_plus_rep1_501bp.bw",
  "GSM10129175_Rec12_WT_plus_rep2_501bp.bw",
  "GSM10129177_Rec12_WT_minus_rep1_501bp.bw",
  "GSM10129178_Rec12_WT_minus_rep2_501bp.bw"
)
chr_name <- c("I","II","III","AB325691","MT","MTR")

#load data
coverage_ls <- lapply(sample_name, function(n){
  
  current_data <- import(n) %>% 
    as.data.table() %>% 
    dplyr::select(seqnames, start, score)
  
  current_data$seqnames <- as.character(current_data$seqnames)
  current_data$start    <- as.numeric(current_data$start)
  current_data$score    <- as.numeric(current_data$score)
  names(current_data) <- c("chr", "coordinate", "count")
  
  return(current_data)
})

#average
replicate_pair  <- list(c(1,2), c(3,4))
ave_sample_name <- c("Rec12_WT_Glucose_plus", "Rec12_WT_Glucose_minus")
ave_coverage_ls <- lapply(replicate_pair, function(n){
  
  ave_count <- lapply(n, function(i){
    coverage_ls[[i]]$count
  }) %>% 
    do.call(cbind, .) %>% 
    rowMeans(na.rm = TRUE)
  
  result <- coverage_ls[[n[1]]]
  result$count <- ave_count
  
  return(result)
})
names(ave_coverage_ls) <- ave_sample_name

#gene data
#pombe_CDS_protein_coding.bed is prepared using BEDtools from gff3 file in pombase.
gene_data <- fread("pombe_CDS_protein_coding.bed")

#plot region
plot_gene_ls <- data.table(name       = c("agl1",  "agl1_2", "fbp1", "inv1"),
                           chr        = c("I",     "I",      "II",   "III"),
                           plot_start = c(2800000, 2950000,  194000, 1719000),
                           plot_end   = c(3100000, 2990000,  204000, 1727000),
                           Rec12_lim  = c(25,      4,        2,      15)
                           )


plot_ls <- lapply(1:nrow(plot_gene_ls), function(n){
  WT_Rec12 <- profile_plot(count_data = ave_coverage_ls, 
                           count_data_col = "count", 
                           chromosome = plot_gene_ls$chr[n], 
                           sample_name_ls = names(ave_coverage_ls), 
                           xlim = c(plot_gene_ls$plot_start[n], plot_gene_ls$plot_end[n]),
                           ylim = c(0, plot_gene_ls$Rec12_lim[n]), 
                           x_label = NULL, 
                           y_label = "Rec12 IP/INPUT", 
                           binned = FALSE, 
                           na_approx = TRUE)
  
  gene <- gene_plot(gene = gene_data,
                    chr = plot_gene_ls$chr[n],
                    xlim = c(plot_gene_ls$plot_start[n], plot_gene_ls$plot_end[n]))
  
  final <- cowplot::align_plots(
    WT_Rec12, gene, 
    align = "v",
    axis = "l"
  )
  
  Fig <- plot_grid(
    plotlist = final,
    ncol = 1
  )
  return(Fig)
})

plot_ls[[1]]
plot_ls[[2]]
plot_ls[[3]]
plot_ls[[4]]
