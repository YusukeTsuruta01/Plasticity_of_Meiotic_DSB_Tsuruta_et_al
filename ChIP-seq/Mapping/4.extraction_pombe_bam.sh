#!/bin/zsh

set -e
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate envs_bowtie2

parameter_dr="$(cd "$(dirname "$0")"; pwd)"
source "${parameter_dr}/0.parameter.sh"

#######################

path_to_data="${base_dr}/bam_file/data/"
output_dr="${base_dr}/bam_file/"

#######################

mkdir -p "${output_dr}/bam_pombe"
mkdir -p "${output_dr}/bam_pombe/data"

for i in "${sample_name[@]}"
do
    echo "Processing sample: ${i}"

    samtools view -b -h -@ 6 "${path_to_data}/${i}_mapped.bam" \
    I II III AB325691 MT MTR | \
    samtools sort -@ 6 -o "${output_dr}/bam_pombe/data/${i}_mapped_pombe.bam"

    samtools index "${output_dr}/bam_pombe/data/${i}_mapped_pombe.bam"

    echo "Completed processing for sample: ${i}"
done
