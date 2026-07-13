#!/usr/bin/env bash
set -euo pipefail

peak_count_matrix="$1"
sample_id="$2"

output_file="${sample_id}_peaks_countMatrix.mm10.bed"

tail -n +3 "$peak_count_matrix" \
| awk -F $'\t' 'BEGIN {OFS = FS} { print $2, $3, $4, $7, $6 }' \
> "$output_file"
