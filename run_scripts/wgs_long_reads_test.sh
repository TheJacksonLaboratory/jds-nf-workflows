#!/bin/bash
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --job-name=wgs_long_reads_mouse
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 36:00:00
#SBATCH --mem=10G
#SBATCH --ntasks=1

cd $SLURM_SUBMIT_DIR

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow/24.10.6

# RUN PIPELINE
nextflow ../main.nf \
--workflow wgs_long_read \
-profile sumner2 \
--gen_org mouse \
--genome_build 'GRCm39' \
--pubdir "/flashscratch/widmas/long_read_outputDir" \
--ref_fa "/projects/compsci/vmp/USERS/widmas/1k_wild_mouse_wgs/metadata/GRCm39_Y_PAR_masked/sequence/GRCm39_masked.fa" \
--ref_fa_indices "/projects/compsci/vmp/USERS/widmas/1k_wild_mouse_wgs/metadata/GRCm39_Y_PAR_masked/indices/GRCm39_masked.fa" \
--minimap2_index "/projects/compsci/vmp/USERS/widmas/1k_wild_mouse_wgs/metadata/GRCm39_Y_PAR_masked/indices/GRCm39_masked.mmi" \
-w "/flashscratch/widmas/long_read_outputDir/work" \
--csv_input "/projects/compsci/vmp/USERS/widmas/1k_wild_mouse_wgs/metadata/simulated_reads_manifest.csv" \
--merge_inds true \
--run_gvcf true \
--comment "This script will run whole genome sequencing analysis on mouse samples using default mm10" \
-resume
