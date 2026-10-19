#!/bin/zsh

set -e
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate envs_bowtie2

parameter_dr="$(cd "$(dirname "$0")"; pwd)"
source "${parameter_dr}/0.parameter.sh"

#######################

path_to_data="${base_dr}/bam_file/bam_SK1/data/"
output_dr="${base_dr}/bam_file/bam_SK1/SK1_peaks/"
SK1_peak_data="${base_dr}/bam_file/bam_SK1/SK1_peaks/merged_peak/merged_peak_SK1.bed"

#######################

mkdir -p "${output_dr}/SK1_peak_read"
mkdir -p "${output_dr}/SK1_peak_read/bam_file"
mkdir -p "${output_dr}/SK1_peak_read/mapping_status"

#extract the mapped read on SK1 hotspot
for i in "${sample_name[@]}"
do
	echo "Processing sample: ${i}"

	bedtools intersect \
    -a "${path_to_data}/${i}_mapped_SK1.bam" \
    -b "${SK1_peak_data}" \
    -ubam | \
    samtools sort -@ 2 -o "${output_dr}/SK1_peak_read/bam_file/${i}_mapped_SK1_peaks.bam"
    samtools index "${output_dr}/SK1_peak_read/bam_file/${i}_mapped_SK1_peaks.bam"

    samtools idxstats "${output_dr}/SK1_peak_read/bam_file/${i}_mapped_SK1_peaks.bam" \
    > "${output_dr}/SK1_peak_read/mapping_status/${i}_mapped_status_SK1_peaks.txt"

    echo "Completed processing for sample: ${i}"
done