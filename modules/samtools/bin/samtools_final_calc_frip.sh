#!/usr/bin/env bash
set -euo pipefail

processed_bam="$1"
reads_peaks_bam="$2"
sample_id="$3"

total_reads=$(samtools view -c "$processed_bam")
reads_in_peaks=$(samtools view -c "$reads_peaks_bam")

frip=$(awk -v rip="$reads_in_peaks" -v total="$total_reads" 'BEGIN { if (total == 0) print 0; else print rip/total }')

echo -e "SAMPLEID\tFRiP\tFiltered Reads\n${sample_id}\t${frip}\t${total_reads}" \
> "${sample_id}_Fraction_reads_in_peak.txt"
