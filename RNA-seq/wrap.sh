#!/bin/zsh

set -e
cd "$(dirname "$0")"

zsh 1.adapter_trimming.sh
zsh 2.mapping_by_hisat2.sh
zsh 3.mapping_status.sh

