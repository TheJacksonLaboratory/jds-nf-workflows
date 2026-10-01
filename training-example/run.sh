#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# On Sumner, uncomment and adjust these lines for the installed Nextflow module.
# module use --append /projects/omics_share/meta/modules
# module load nextflow

nextflow run main.nf \
    -profile local \
    --csv_input "${SCRIPT_DIR}/examples/messages.csv" \
    --pubdir "${SCRIPT_DIR}/results"

cat "${SCRIPT_DIR}/results"/*.summary.txt
cat "${SCRIPT_DIR}/results/all_summaries.txt"