#!/bin/zsh

set -e
cd "$(dirname "$0")"

zsh 1.adapter_trimming.sh
zsh 2.mapping_by_bowtie2.sh
zsh 3.mapping_status.sh
zsh 4.extraction_pombe_bam.sh
zsh 5.extraction_SK1_bam.sh
zsh 6.peak_call_for_SK1.sh
zsh 7.merge_SK1_peak.sh
zsh 8.extract_mapped_read_SK1_peak.sh
zsh 9.coverage_from_bam_file.sh

