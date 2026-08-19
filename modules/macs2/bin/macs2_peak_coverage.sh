#!/usr/bin/env bash
set -euo pipefail

narrow_peaks="$1"
sample_id="$2"
output_file="${sample_id}_peaks.narrowPeak.saf"

awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' "$narrow_peaks" \
> "$output_file"
