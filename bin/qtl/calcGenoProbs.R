#!/usr/bin/env Rscript
library(qtl2)
library(parallel)

################################################################################
# Perform initial calculations of genotype probabilities.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20241216
################################################################################ 

args <- commandArgs(trailingOnly = TRUE)
# 
# test_dir <- "/flashscratch/widmas/HR_QC_outputDir/work/b6/8d72b12fb2ba334a2a162647f22aa9"
# setwd(test_dir)

# cross object from WRITE_CROSS
# cross_file <- args[1]

# Load in DO cross data
cross <- readRDS(list.files(pattern = "preQC"))

# Calculate genotype probs
pr <- qtl2::calc_genoprob(cross = cross, 
                          map = cross$pmap, 
                          error_prob = 0.002, 
                          lowmem = FALSE,
                          cores = (parallel::detectCores()/2), 
                          quiet = F)
saveRDS(pr, file = "pr_36state.rds")