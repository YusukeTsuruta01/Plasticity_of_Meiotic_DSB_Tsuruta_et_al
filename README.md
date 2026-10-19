# Rec12 ChIP-seq analysis, Tsuruta *et al*.

This repository contains the shell and R scripts for the analysis in Tsuruta *et al*. (doi: xxxx)\
from Kouji Hirota laboratory in Tokyo Metropolitan University.\
The details of analysis are written in the method section of paper.\
The data used for analysis are available at GEO and supplementary table of paper (Table_S3_NGS_data_in_hotspots.xlsx).\
Further information and requests for resources should be directed to and will be fulfilled by the lead contact, Kouji Hirota (khirota@tmu.ac.jp).

## ChIP-seq
Raw data and bigwig file of ChIP-seq are available at GEO (GSE351198).
### Mapping
0.parameter.sh\
1.adapter_trimming.sh\
2.mapping_by_bowtie2.sh\
3.mapping_status.sh\
4.extraction_pombe_bam.sh\
5.extraction_SK1_bam.sh\
6.peak_call_for_SK1.sh\
7.merge_SK1_peak.sh\
8.extract_mapped_read_SK1_peak.sh\
9.coverage_from_bam_file.sh\
wrap.sh
### Normalization, calculation of IP/INPUT, and smoothing
10.binning_10bp.R\
11.mask_based_on_INPUT.R\
12.normalize_and_IP_vs_INPUT.R\
13.smoothing.R\
14.average.R
### Peak calling
15.peak_call.R\
16.merge_peaks.R\
17.merge_adjacent_peak_500bp.sh
### Analysis of peak details
18.calculation_of_peak_intensity.R\
19.judge_overlap_with_gene.R
### Other analysis
20.multiple_regression_analysis.R
## RNA-seq
Raw data of RNA-seq and TPM normalized count data are available at GEO (GSE351199).
### Mapping
0.parameter.sh\
1.adapter_trimming.sh\
2.mapping_by_hisat2.sh\
3.mapping_status.sh\
wrap.sh
## Figure generation
1.bar_plot.R\
2.line_plot.R\
3.ChIP_seq_profile_plot.R\
4.scatter_plot.R\
5.violin_plot.R\
6.correlation_matrix.R\
7.MA_plot.R
## Statistical analysis
1.analysis_using_limma.R\
2.ANOVA_and_TukeyHSD.R\
3.Kruskal_Wallis_and_Dunn.R\
4.deference_of_correlation_using_cocor.R
