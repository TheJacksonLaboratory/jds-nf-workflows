#!/bin/bash
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --job-name=wgs_nf-test
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 12:00:00
#SBATCH --mem=50G
#SBATCH --ntasks=1

# nf-test container
NF_TEST_CONTAINER="/projects/compsci/vmp/USERS/widmas/containers/nf-test_nextflow.img"

# LOAD SINGULARITY
module load singularity

# RUN NF-TEST
singularity exec ${NF_TEST_CONTAINER} nf-test test tests/workflows/wgs.nf.test