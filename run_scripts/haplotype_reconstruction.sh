#!/bin/bash
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --job-name=haplotype_reconstruction
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
--workflow haplotype_reconstruction \
-profile sumner2 \
--pubdir "/flashscratch/widmas/hr_outputDir" \
-w "/flashscratch/widmas/hr_outputDir/work" \
--csv_input "/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/sample_sheets/haplotype_reconstruction_input.csv" \
--rerun true \
--correct_ids false \
--remove_markers false \
--comment "This script will run haplotype reconstruction on mouse genotyped using GigaMUGA on default mm10 coordinates"
