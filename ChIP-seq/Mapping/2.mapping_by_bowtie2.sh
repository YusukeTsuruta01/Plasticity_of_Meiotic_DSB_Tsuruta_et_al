#!/bin/zsh

set -e
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate envs_bowtie2

parameter_dr="$(cd "$(dirname "$0")"; pwd)"
source "${parameter_dr}/0.parameter.sh"

#######################

path_to_data="${base_dr}/raw_data/adaptor_trimmed_data/data/"
output_dr="${base_dr}/"

#######################

mkdir -p "${output_dr}/bam_file"
mkdir -p "${output_dr}/bam_file/data"
mkdir -p "${output_dr}/bam_file/mapping_log"

#mapping by bowtie2
for i in "${sample_name[@]}"
do
	echo "Processing sample: ${i}"
	echo "Running bowtie2..."
	
	{ bowtie2 -x "$index_name" \
	-1 "${path_to_data}/${i}_1_trimmed.fastq.gz" \
	-2 "${path_to_data}/${i}_2_trimmed.fastq.gz" \
	-p 6 -N 0 ; }\
	2> "${output_dr}/bam_file/mapping_log/${i}_bowtie2.log" | \
	grep -v "XS:" | \
	samtools view -b -f 2 -F 256 -@ 6 | \
	samtools sort -@ 6 -o "${output_dr}/bam_file/data/${i}_mapped.bam"

	echo "Indexing bam file..."
	samtools index "${output_dr}/bam_file/data/${i}_mapped.bam"

  echo "Completed processing for sample: ${i}"

done
