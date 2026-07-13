#!/usr/bin/env bash
set -euo pipefail

target_bed="$1"
bam_file="$2"
sample_id="$3"
output_file="${sample_id}_amplicon_coverage_metrics.txt"

# Total bases that map/align to the on-target region.
bases_on_target=$(coverageBed -a "$target_bed" -b "$bam_file" | awk '{if($7>0) total+=$7}END{print total}')

# Total length covered by BAM alignment.
total_bases_covered=$(genomeCoverageBed -ibam "$bam_file" -bg | awk '{if($4>0) total += ($3-$2)}END{print total}')

# On-target percent.
awk -v a="$bases_on_target" -v b="$total_bases_covered" 'BEGIN { printf "on_target_percent\t%s\n", (a/b)*100 }' </dev/null > "$output_file"

# Depth of 20% of the average coverage.
perc_mean=$(coverageBed -d -a "$target_bed" -b "$bam_file" | awk '{if($7>0) total+=1;s+=$7}END{print (s/total)*.2}')

# Total capture array bases.
total_target_bases=$(awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' "$target_bed")

# Coverage uniformity.
coverageBed -d -a "$target_bed" -b "$bam_file" | awk -v percmean="$perc_mean" -v totalbases="$total_target_bases" '{if($7>percmean) total+=1;s+=$7}END{print "coverage_uniformity\t"(total/totalbases)*100}' >> "$output_file"
