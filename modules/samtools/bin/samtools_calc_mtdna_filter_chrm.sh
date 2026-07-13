#!/usr/bin/env bash
set -euo pipefail

rmdup_bam_file="$1"
sample_id="$2"
mt_name="$3"
cpus="$4"

# Get mitochondrial read counts from BAM.
mt_reads=$(samtools idxstats "$rmdup_bam_file" | awk -v mt="$mt_name" '$1==mt{print $3}')

# Get total read counts from BAM.
total_reads=$(samtools idxstats "$rmdup_bam_file" | awk '{SUM += $3} END {print SUM}')

if [[ -z "$mt_reads" ]] || [[ "$mt_reads" == "" ]]; then
    mt_reads=0
fi

if [[ -z "$total_reads" ]] || [[ "$total_reads" == "0" ]]; then
    total_reads=1
fi

# Calculate %mtDNA and write report.
echo -e "sampleID\tPerc mtDNA\n${sample_id}\t$(bc <<< \"scale=2;100*${mt_reads}/${total_reads}\")" > "${sample_id}_mtDNA_Content.txt"

# Filter mitochondrial reads from BAM and index output.
samtools view -@ "$cpus" -h "$rmdup_bam_file" \
| grep -v "$mt_name" \
| samtools sort -@ "$cpus" -O bam -o "${sample_id}.sorted.rmDup.rmChrM.bam"

samtools index "${sample_id}.sorted.rmDup.rmChrM.bam"
