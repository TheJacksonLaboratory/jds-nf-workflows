#!/bin/bash
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --job-name=haplotype_reconstruction
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
-profile sumner2 \
--workflow haplotype_reconstruction \
--csv_input  "/projects/compsci/vmp/USERS/widmas/SODO/HR_rerun_input.csv" \
--rerun false \
--correct_ids true \
--remove_markers true \
--pubdir "/flashscratch/widmas/sodo_hr" \
-w "/flashscratch/widmas/sodo_hr/work" \
--comment "This script will run haplotype reconstruction on mouse genotyped using GigaMUGA on default mm10 coordinates"
