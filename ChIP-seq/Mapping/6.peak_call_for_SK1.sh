#!/bin/zsh
set -e
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate envs_macs3

parameter_dr="$(cd "$(dirname "$0")"; pwd)"
source "${parameter_dr}/0.parameter.sh"

#######################

#sample name
INPUT_sample_name=(
	"WT_glu_plus_N1_IN"
	"WT_glu_minus_N1_IN"
	"WT_glu_plus_N2_IN"
	"WT_glu_minus_N2_IN"
)
IP_sample_name=(
	"WT_glu_plus_N1_IP"
	"WT_glu_minus_N1_IP"
	"WT_glu_plus_N2_IP"
	"WT_glu_minus_N2_IP"
)
output_file_name=(
    "WT_glu_plus_N1"
	"WT_glu_minus_N1"
	"WT_glu_plus_N2"
	"WT_glu_minus_N2"
)

path_to_data="${base_dr}/bam_file/bam_SK1/data/"
output_dr="${base_dr}/bam_file/bam_SK1/"
genome_size=12631379

#option
q_val=0.05

#######################

mkdir -p "${output_dr}/SK1_peaks"
mkdir -p "${output_dr}/SK1_peaks/peak_data"

#peak calling by macs3
INPUT_sample_number=${#INPUT_sample_name[@]}

for i in $(seq 1 ${INPUT_sample_number})
do
	echo "Processing sample: ${output_file_name[i]}"
    macs3 callpeak \
    -t "${path_to_data}/${IP_sample_name[i]}_mapped_SK1.bam" \
    -c "${path_to_data}/${INPUT_sample_name[i]}_mapped_SK1.bam" \
    --outdir "${output_dr}/SK1_peaks/peak_data" \
    --nomodel \
    -q ${q_val} \
    -f BAMPE \
    -g ${genome_size} \
    -n "${output_file_name[i]}"
	
	echo "Completed processing for sample: ${output_file_name[i]}"
done
