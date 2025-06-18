#!/usr/bin/env Rscript

################################################################################
# Execute permutation tests for QTL mapping
# - Perform permutation tests for each phenotype, using the requested covariates.
# - Write out permutation results.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20250618
################################################################################

library(qtl2)
args <- commandArgs(trailingOnly = TRUE)
covar_file        <- args[1]
cross_file        <- args[2]
genoprobs_file    <- args[3]
alleleprobs_file  <- args[4]
kinship_file      <- args[5]
pheno_file        <- args[6]
n_perms           <- as.numeric(args[7])

# Test files
# testDir <- "/flashscratch/widmas/qtl_mapping_qc"
# alleleprobs_file  <- file.path(testDir,"apr.rds")
# cross_file        <- file.path(testDir,"cross.rds")
# kinship_file      <- file.path(testDir,"kinship.rds")
# pheno_file        <- file.path(testDir,"Ins_tAUC_pheno.csv")
# n_perms <- 5

# Read in files
alleleprobs <- readRDS(alleleprobs_file)
cross <- readRDS(cross_file)
kinship <- readRDS(kinship_file)
pheno <- read.csv(pheno_file, row.names = 1)

# make covar matrix
addcovar = model.matrix(~ sex + gen, data = cross$covar)[,-1,drop = FALSE]

# Run permutation tests
perms <- qtl2::scan1perm(genoprobs = alleleprobs,
                         pheno = pheno, 
                         kinship = kinship,
                         addcovar = addcovar,
                         n_perm = n_perms,
                         perm_Xsp=TRUE,
                         chr_lengths=chr_lengths(cross$pmap),
                         cores = parallel::detectCores())

# Save the files
saveRDS(perms, file = paste0(colnames(pheno),"_scan1perms.rds"))
