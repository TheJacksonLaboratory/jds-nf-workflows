#!/usr/bin/env Rscript

################################################################################
# Run QTL mapping
# - Map each phenotype, using the requested covariates.
# - Write out LODs and make QTL plots.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20250716
################################################################################

library(qtl2)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
covar_file        <- args[1]
cross_file        <- args[2]
genoprobs_file    <- args[3]
alleleprobs_file  <- args[4]
kinship_file      <- args[5]
pheno_file        <- args[6]
covar_info_file   <- args[7]

# # Test files
# testDir <- "/flashscratch/widmas/qtl_mapping_outputDir/work/4d/29dabfc985762e9992360a89f26f2b"
# setwd(testDir)
# alleleprobs_file  <- "apr.rds"
# cross_file        <- list.files(testDir, pattern = "cross.rds")
# kinship_file      <- "kinship.rds"
# pheno_file        <- list.files(testDir, pattern = "pheno.csv")
# covar_info_file   <- list.files(testDir, pattern = "covar_info.csv")


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
  
  # QTL scan with interaction
  scan1out <- qtl2::scan1(genoprobs = alleleprobs,
                          pheno = pheno, 
                          kinship = kinship,
                          addcovar = covar_matrix,
                          intcovar = interactive_covariate,
                          cores = parallel::detectCores())
  
} else {
  
  # QTL scan without interaction
  scan1out <- qtl2::scan1(genoprobs = alleleprobs,
                          pheno = pheno, 
                          kinship = kinship,
                          addcovar = covar_matrix,
                          cores = parallel::detectCores())
  
}


# Plot LODs
png(paste0(colnames(pheno),"_scan1.png"))
qtl2::plot_scan1(x = scan1out, map = cross$pmap, main = colnames(pheno))
dev.off()

# Save the files
saveRDS(scan1out, file = paste0(colnames(pheno),"_scan1out.rds"))
saveRDS(cross, file = paste0(colnames(pheno),"_cross.rds"))
