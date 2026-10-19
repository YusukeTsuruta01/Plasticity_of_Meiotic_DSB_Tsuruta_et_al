#!/bin/zsh

set -e
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate envs_bowtie2

parameter_dr="$(cd "$(dirname "$0")"; pwd)"
source "${parameter_dr}/0.parameter.sh"

#######################

path_to_data="${base_dr}/bam_file/bam_pombe/data/"
output_dr="${base_dr}/"

#######################

mkdir -p "${output_dr}/coverage_data"
mkdir -p "${output_dr}/coverage_data/from_bam_file"

#Get coverage data from bam file
for i in "${sample_name[@]}"
do
	echo "Processing sample: ${i}"
	samtools depth -aa "${path_to_data}/${i}_mapped_pombe.bam" | \
	gzip \
	> "${output_dr}/coverage_data/from_bam_file/${i}_coverage_map.txt.gz"

	echo "Completed processing for sample: ${i}"
done






