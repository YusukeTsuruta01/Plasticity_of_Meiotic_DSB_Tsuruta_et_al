#!/bin/zsh

set -e
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate envs_bowtie2

parameter_dr="$(cd "$(dirname "$0")"; pwd)"
source "${parameter_dr}/0.parameter.sh"

#######################

path_to_data="${base_dr}/raw_data/data/"
output_dr="${base_dr}/raw_data/"

#######################

mkdir -p "${output_dr}/adaptor_trimmed_data"
mkdir -p "${output_dr}/adaptor_trimmed_data/data"
mkdir -p "${output_dr}/adaptor_trimmed_data/trimming_report"
mkdir -p "${output_dr}/adaptor_trimmed_data/QC"

#adaptor trimming by fastp
sample_number=${#original_name[@]}
echo "Total sample number: ${sample_number}"

for i in $(seq 1 ${sample_number})
do
	echo "Processing sample: ${original_name[i]}"
	
	fastp \
	-i "${path_to_data}/${original_name[i]}_1.fq.gz" \
	-I "${path_to_data}/${original_name[i]}_2.fq.gz" \
	-o "${output_dr}/adaptor_trimmed_data/data/${sample_name[i]}_1_trimmed.fastq.gz" \
	-O "${output_dr}/adaptor_trimmed_data/data/${sample_name[i]}_2_trimmed.fastq.gz" \
	-q 20 -l 20 -3 \
	-h "${output_dr}/adaptor_trimmed_data/trimming_report/${sample_name[i]}_trimming_report.html" \
	-j "${output_dr}/adaptor_trimmed_data/trimming_report/${sample_name[i]}_trimming_report.json"
	
	echo "Completed processing for sample: ${original_name[i]}"
done

#check by fastqc
fastqc "${output_dr}"/adaptor_trimmed_data/data/*.fastq.gz -o "${output_dr}/adaptor_trimmed_data/QC" -t 6
