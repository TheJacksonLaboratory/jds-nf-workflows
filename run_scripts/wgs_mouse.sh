#!/bin/bash
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --job-name=wgs_mouse_WMGP
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH --mem=5G
#SBATCH --ntasks=1

cd $SLURM_SUBMIT_DIR

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow/24.10.6

# RUN PIPELINE
nextflow ../main.nf \
--workflow wgs \
-profile sumner2 \
--gen_org mouse \
--genome_build 'GRCm39' \
--ref_fa "/projects/compsci/vmp/USERS/widmas/1k_wild_mouse_wgs/metadata/GRCm39_Y_PAR_masked/sequence/GRCm39_masked.fa" \
--ref_fa_indices "/projects/compsci/vmp/USERS/widmas/1k_wild_mouse_wgs/metadata/GRCm39_Y_PAR_masked/sequence/GRCm39_masked.fa.fai" \
--ref_fa_dict "/projects/compsci/vmp/USERS/widmas/1k_wild_mouse_wgs/metadata/GRCm39_Y_PAR_masked/combined_ref_set/GRCm39_masked.dict" \
--csv_input "/projects/compsci/vmp/USERS/widmas/1k_wild_mouse_wgs/metadata/20251211_wgs_from_bam_manifest.csv" \
--bam_input true \
--deepvariant false \
--run_sv true \
--run_gvcf true \
--pubdir "/flashscratch/widmas/wmgp_gatk_sv" \
-w "/flashscratch/widmas/wmgp_gatk_sv_mito/work" \
--comment "This script will run whole genome sequencing analysis on mouse samples using default mm10" \
-resume
