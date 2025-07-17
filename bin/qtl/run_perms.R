#!/usr/bin/env Rscript

################################################################################
# Execute permutation tests for QTL mapping
# - Perform permutation tests for each phenotype, using the requested covariates.
# - Write out permutation results.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20250716
################################################################################

library(qtl2)
args <- commandArgs(trailingOnly = TRUE)
covar_file        <- args[1]
cross_file        <- args[2]
genoprobs_file    <- args[3]
alleleprobs_file  <- args[4]
kinship_file      <- args[5]
pheno_file        <- args[6]
covar_info_file   <- args[7]
n_perms           <- as.numeric(args[8])

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
covar_info <- read.csv(covar_info_file)

# end if the phenotypes listed in different files don't match
stopifnot(colnames(pheno) == unique(covar_info$phenotype))

# make covar matrix
covar_formula <- paste0(covar_info$covar, collapse="+")
covar_formula <- paste0("~", covar_formula)
covar_matrix <- stats::model.matrix.lm(
  stats::as.formula(covar_formula),
  data = cross$covar,
  na.action = stats::na.pass
)

# drop the intercept column
covar_matrix <- covar_matrix[, -1, drop = FALSE]

# drop the covar column if it has all identical values
covar_matrix <- covar_matrix[, apply(covar_matrix, 2, function(col) length(unique(col)) > 1), drop = FALSE]

# detect any interactive covariates
if(any(covar_info$interactive)){
  intcovar <- covar_info$covar[which(covar_info$interactive)]
  interactive_covariate <-
    covar_matrix[, which(grepl(intcovar, colnames(covar_matrix), ignore.case = T))]
  
  # Run permutation tests without interaction
  perms <- qtl2::scan1perm(genoprobs = alleleprobs,
                           pheno = pheno, 
                           kinship = kinship,
                           addcovar = covar_matrix, 
                           intcovar = interactive_covariate,
                           n_perm = n_perms,
                           perm_Xsp=TRUE,
                           chr_lengths=chr_lengths(cross$pmap),
                           cores = parallel::detectCores())
  
} else {
  
  # Run permutation tests without interaction
  perms <- qtl2::scan1perm(genoprobs = alleleprobs,
                           pheno = pheno, 
                           kinship = kinship,
                           addcovar = covar_matrix,
                           n_perm = n_perms,
                           perm_Xsp=TRUE,
                           chr_lengths=chr_lengths(cross$pmap),
                           cores = parallel::detectCores())
  
}


# Save the files
saveRDS(perms, file = paste0(colnames(pheno),"_scan1perms.rds"))
