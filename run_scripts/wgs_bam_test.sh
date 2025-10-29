#!/bin/bash
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --job-name=wgs_mouse
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH --mem=100G
#SBATCH --ntasks=1

cd $SLURM_SUBMIT_DIR

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow/24.10.6

# RUN PIPELINE
nextflow ../main.nf \
--workflow wgs \
-profile sumner2 \
--gen_org human \
--genome_build 'GRCh38' \
--pubdir "/flashscratch/widmas/wgsBam_human_outputDir" \
-w "/flashscratch/widmas/wgsBam_human_outputDir/work" \
--csv_input "/projects/compsci/vmp/USERS/widmas/jds-nf-test/wgs/human/hg38_WGS_input.csv" \
--comment "This script will run whole genome sequencing analysis on mouse samples using default mm10" \
-resume
