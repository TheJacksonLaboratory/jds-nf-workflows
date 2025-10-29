#!/bin/bash
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --job-name=qtl_mapping
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
--workflow qtl_mapping \
-profile sumner2 \
--pubdir "/flashscratch/widmas/qtl_mapping_outputDir" \
-w "/flashscratch/widmas/qtl_mapping_outputDir/work" \
--csv_input "" \
--n_perms 5 \
--comment "This script will run QTL mapping on mouse samples using default mm10 coordinates" \
-resume
