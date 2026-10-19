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

mkdir -p "${output_dr}/bam_SK1"
mkdir -p "${output_dr}/bam_SK1/data"

for i in "${sample_name[@]}"
do
    echo "Processing sample: ${i}"

    samtools view -b -h -@ 6 "${path_to_data}/${i}_mapped.bam" \
    SK1_1 SK1_2 SK1_3 SK1_4 SK1_5 SK1_6 SK1_7 SK1_8 SK1_9 SK1_10 SK1_11 SK1_12 SK1_13 SK1_14 SK1_15 SK1_16 SK1_MT | \
    samtools sort -@ 6 -o "${output_dr}/bam_SK1/data/${i}_mapped_SK1.bam"

    samtools index "${output_dr}/bam_SK1/data/${i}_mapped_SK1.bam"

    echo "Completed processing for sample: ${i}"
done
