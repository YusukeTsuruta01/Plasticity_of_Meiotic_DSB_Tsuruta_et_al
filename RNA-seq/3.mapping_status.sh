#!/bin/zsh

set -e
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate envs_hisat2

parameter_dr="$(cd "$(dirname "$0")"; pwd)"
source "${parameter_dr}/0.parameter.sh"

#######################

path_to_data="${base_dr}/bam_file/data/"
output_dr="${base_dr}/"

#######################

mkdir -p "${path_to_data}/mapping_status"

#check mapping status
for i in "${sample_name[@]}"
do
    samtools idxstats "${path_to_data}/${i}_mapped.bam" > "${path_to_data}/mapping_status/${i}_mapped_status.txt"
done

