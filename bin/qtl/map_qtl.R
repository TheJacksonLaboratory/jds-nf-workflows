#!/usr/bin/env Rscript

################################################################################
# Run QTL mapping
# - Map each phenotype, using the requested covariates.
# - Write out LODs and make QTL plots.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20260121
################################################################################

library(qtl2)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
covar_file        <- args[1]
map_file          <- args[2]
genoprobs_file    <- args[3]
alleleprobs_file  <- args[4]
kinship_file      <- args[5]
pheno_file        <- args[6]
covar_info_file   <- args[7]
n_cores           <- as.numeric(args[8])

# Read in files
alleleprobs <- readRDS(alleleprobs_file)
map <- readRDS(map_file)
kinship <- readRDS(kinship_file)
pheno <- read.csv(pheno_file, row.names = 1)
covar_info <- read.csv(covar_info_file)
covar <- read.csv(covar_file)

# end if the phenotypes listed in different files don't match
stopifnot(colnames(pheno) == unique(covar_info$phenotype))

# make covar matrix
covar_formula <- paste0(covar_info$covar, collapse="+")
covar_formula <- paste0("~", covar_formula)
covar_matrix <- stats::model.matrix.lm(
  stats::as.formula(covar_formula),
  data = covar,
  na.action = stats::na.pass
)

# drop the intercept column
covar_matrix <- covar_matrix[, -1, drop = FALSE]

# drop the covar column if it has all identical values
covar_matrix <- covar_matrix[, apply(covar_matrix, 2, function(col) length(unique(col)) > 1), drop = FALSE]
rownames(covar_matrix) <- covar$id

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
                          cores = n_cores)
  
} else {
  
  # QTL scan without interaction
  scan1out <- qtl2::scan1(genoprobs = alleleprobs,
                          pheno = pheno, 
                          kinship = kinship,
                          addcovar = covar_matrix,
                          cores = n_cores)
  
}


# Plot LODs
png(paste0(colnames(pheno),"_scan1.png"))
qtl2::plot_scan1(x = scan1out, map = map, main = colnames(pheno))
dev.off()

# Save the files
saveRDS(scan1out, file = paste0(colnames(pheno),"_scan1out.rds"))
saveRDS(map, file = paste0(colnames(pheno),"_map.rds"))
