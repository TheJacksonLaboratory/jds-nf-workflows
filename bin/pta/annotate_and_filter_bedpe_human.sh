#!/bin/bash

# Usage:
# ./annotate_and_filter_bedpe.sh <input_bedpe> <input_cnv> <max_changepoint_distance> <filter_databases> <genome> <out_somatic> <out_highconf>

# Arguments
BEDPE_IN="$1"
CNV_IN="$2"
OUT_SOMATIC="$4"
OUT_HIGHCONF="$5"

# Constants
MAX_CP_DIST="1000"
FILTER_DBS="DGV,1000G,PON"
GENOME="GRCh38"
ENSEMBL_BED="/projects/omics_share/human/GRCh38/supporting_files/PTA_inputs/annotations/ensembl_genes_unique_sorted.final.v93.chr.sorted.bed"         # Update with actual path or pass as argument
CANCER_CENSUS_BED="/projects/omics_share/human/GRCh38/supporting_files/PTA_inputs/annotations/cancer_gene_census.GRCh38-v92.bed"    # Update with actual path or pass as argument

# Apptainer image
APPTAINER_IMG="/projects/compsci/omics_share/meta/containers/quay.io-jaxcompsci-r-sv_cnv_annotate-4.1.1.img"

# Step 1: Annotate BEDPE with gene information
GENE_ANNOTATE_R="/projects/omics_share/meta/benchmarking/ngs-ops-nf-pipelines/bin/pta/annotate-bedpe-with-genes.r"
GENE_ANNOTATED_BEDPE="gene_annotated_${BEDPE_IN##*/}"

apptainer run "$APPTAINER_IMG" Rscript "$GENE_ANNOTATE_R" \
    --ensembl="$ENSEMBL_BED" \
    --cancer_census="$CANCER_CENSUS_BED" \
    --bedpe="$BEDPE_IN" \
    --out_file="$GENE_ANNOTATED_BEDPE"

# Step 2: Annotate BEDPE with CNV changepoints
ANNOTATE_R="/projects/omics_share/meta/benchmarking/ngs-ops-nf-pipelines/bin/pta/annotate-bedpe-with-cnv.r"
CNV_ANNOTATED_BEDPE="cnv_gene_annotated_${BEDPE_IN##*/}"

apptainer run "$APPTAINER_IMG" Rscript "$ANNOTATE_R" \
  --bedpe "$GENE_ANNOTATED_BEDPE" \
  --cnv "$CNV_IN" \
  --out_file "$CNV_ANNOTATED_BEDPE"

# Step 3: Filter annotated BEDPE for somatic and high-confidence variants
FILTER_R="/projects/omics_share/meta/benchmarking/ngs-ops-nf-pipelines/bin/pta/filter-bedpe.r"

apptainer run "$APPTAINER_IMG" Rscript "$FILTER_R" \
  --bedpe "$CNV_ANNOTATED_BEDPE" \
  --max_changepoint_distance "$MAX_CP_DIST" \
  --filter_databases "$FILTER_DBS" \
  --out_file_somatic "$OUT_SOMATIC" \
  --out_file_highconf "$OUT_HIGHCONF" \
  --genome "$GENOME"
