#!/bin/zsh
set -e
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate envs_bedtools

parameter_dr="$(cd "$(dirname "$0")"; pwd)"
source "${parameter_dr}/0.parameter.sh"

#######################

path_to_data="${base_dr}/bam_file/bam_SK1/SK1_peaks/peak_data"
output_dr="${base_dr}/bam_file/bam_SK1/SK1_peaks"
output_name="merged_peak_SK1.bed"

#######################

#merge SK1 peaks
mkdir -p "${output_dr}/merged_peak"

awk '
    BEGIN {
        OFS = "\t"
        order["SK1_1"] = 1;
        order["SK1_2"] = 2;
        order["SK1_3"] = 3;
        order["SK1_4"] = 4;
        order["SK1_5"] = 5;
        order["SK1_6"] = 6;
        order["SK1_7"] = 7;
        order["SK1_8"] = 8;
        order["SK1_9"] = 9;
        order["SK1_10"] = 10;
        order["SK1_11"] = 11;
        order["SK1_12"] = 12;
        order["SK1_13"] = 13;
        order["SK1_14"] = 14;
        order["SK1_15"] = 15;
        order["SK1_16"] = 16;
        order["SK1_MT"] = 17
    }
    { print order[$1], $0 } 
    ' ${path_to_data}/*.narrowPeak | \
sort -k1,1n -k3,3n | \
cut -f2- | \
bedtools merge | \
awk 'BEGIN {OFS = "\t"}{ print $1, $2, $3, "SK1_peaks_" NR}' \
> "${output_dr}/merged_peak/${output_name}"

