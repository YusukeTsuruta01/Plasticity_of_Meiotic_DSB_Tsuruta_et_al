#!/bin/zsh
set -e
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate envs_bedtools

base_dr=""

path_to_data="${base_dr}/analysis/peak_call/merged_peak/"
output_dr="${base_dr}/analysis/peak_call/merged_peak/"

sample_name=(
    "WT_glu_plus_minus_rec12_merged.bed"
)

peak_name=(
    "rec12_peaks"
)

output_name=(
    "rec12_peaks_500bp_merged"
)

merge_threshold=500
col_of_peak_name=4

mkdir -p "${output_dr}/merge_adjacent_${merge_threshold}bp"
#merge peaks in close ploximity in each sample
sample_number=${#sample_name[@]}
for i in $(seq 1 ${sample_number})
do
    awk '
    BEGIN {
        OFS = "\t"
        order["I"]   = 1;
        order["II"]  = 2;
        order["III"] = 3;
    }
    { print order[$1], $0 } 
    ' "${path_to_data}/${sample_name[i]}" | \
    sort -k1,1n -k3,3n | \
    cut -f2- | \
    bedtools merge \
    -d ${merge_threshold} \
    -c ${col_of_peak_name} \
    -o collapse |\
    awk -v current_peak_name="${peak_name[i]}_" \
    'BEGIN {OFS = "\t"}{ print $1, $2, $3, current_peak_name NR, $4 }' \
    > "${output_dr}/merge_adjacent_${merge_threshold}bp/${output_name[i]}_merge_adjacent_${merge_threshold}bp.bed"
done

