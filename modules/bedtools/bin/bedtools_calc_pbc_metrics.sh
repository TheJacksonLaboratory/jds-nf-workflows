#!/usr/bin/env bash
set -euo pipefail

tmp_bam="$1"
sample_id="$2"
output_file="${sample_id}.pbc.qc"

if bedtools bamtobed -bedpe -i "$tmp_bam" \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' \
    | grep -v 'MT' \
    | sort \
    | uniq -c \
    | awk -v sample="$sample_id" 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "SAMPLEID\tMT\tM0\tM1\tM2\tNRF\tPBC1\tPBC2\n%s\t%d\t%d\t%d\t%d\t%f\t%f\t%f\n",sample,mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' \
    > "$output_file"; then
    :
else
    echo -e "NA\tNA\tNA\tNA\tNA\tNA\tNA" > "$output_file"
fi
