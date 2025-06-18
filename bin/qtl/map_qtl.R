#!/usr/bin/env Rscript

################################################################################
# Run QTL mapping
# - Map each phenotype, using the requested covariates.
# - Write out LODs and make QTL plots.
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

# Test files
# testDir <- "/flashscratch/widmas/qtl_mapping_qc"
# alleleprobs_file  <- file.path(testDir,"apr.rds")
# cross_file        <- file.path(testDir,"cross.rds")
# kinship_file      <- file.path(testDir,"kinship.rds")
# pheno_file        <- file.path(testDir,"Ins_tAUC_pheno.csv")


# Read in files
alleleprobs <- readRDS(alleleprobs_file)
cross <- readRDS(cross_file)
kinship <- readRDS(kinship_file)
pheno <- read.csv(pheno_file, row.names = 1)

# make covar matrix
addcovar = model.matrix(~ sex + gen, data = cross$covar)[,-1,drop = FALSE]

# QTL scan
scan1out <- qtl2::scan1(genoprobs = alleleprobs,
                        pheno = pheno, 
                        kinship = kinship,
                        addcovar = addcovar, 
                        cores = parallel::detectCores())

# Plot LODs
png(paste0(colnames(pheno),"_scan1.png"))
qtl2::plot_scan1(x = scan1out, map = cross$pmap)
dev.off()

# Save the files
saveRDS(scan1out, file = paste0(colnames(pheno),"_scan1out.rds"))
