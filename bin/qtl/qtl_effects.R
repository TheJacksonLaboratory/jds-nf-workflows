#!/usr/bin/env Rscript

################################################################################
# Estimate QTL effects
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20251031
################################################################################

library(qtl2)
library(dplyr)

# Test files
args <- commandArgs(trailingOnly = TRUE)
phenotype <- args[1]
alleleprobs_file <- args[2]
kinship_file <- args[3]
covar_file <- args[4]
phenotype_file <- args[5]
covar_info_file <- args[6]
map_file <- args[7]
chrom <- args[8]
pos_peak <- as.numeric(args[9])
pos_start <- as.numeric(args[10])
pos_end <- as.numeric(args[11])
scan1out_file <- args[12]

# Read in the files
alleleprobs <- readRDS(alleleprobs_file)
kinship <- readRDS(kinship_file)
covar <- read.csv(covar_file, row.names = 1)
pheno <- read.csv(phenotype_file, row.names = 1)
covar_info <- read.csv(covar_info_file)
map <- readRDS(map_file)
scan1out <- readRDS(scan1out_file)

# Subset probs to QTL location
qtl_probs <- qtl2::pull_genoprobint(genoprobs = alleleprobs, 
                                    map = map, 
                                    chr = chrom, 
                                    interval = c(pos_start, pos_end))
qtl_peak_probs <- qtl2::pull_genoprobpos(genoprobs = alleleprobs, 
                                         map = map, 
                                         chr = chrom, 
                                         pos = pos_peak)

# Make covar matrix
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

# detect any interactive covariates
if(any(covar_info$interactive)){
  intcovar <- covar_info$covar[which(covar_info$interactive)]
  interactive_covariate <-
    covar_matrix[, which(grepl(intcovar, colnames(covar_matrix), ignore.case = T))]
  
  # Pull effects
  scan1coef_out <- qtl2::scan1coef(genoprobs = qtl_probs, 
                                   pheno = pheno, 
                                   kinship = kinship[[chrom]], 
                                   addcovar = covar_matrix, 
                                   intcovar = interactive_covariate)
  # Pull blups
  scan1blup_out <- qtl2::scan1blup(genoprobs = qtl_probs, 
                                   pheno = pheno, 
                                   kinship = kinship[[chrom]], 
                                   addcovar = covar_matrix)
  
  # Fit effects at peak
  fit1_out <- qtl2::fit1(genoprobs = qtl_peak_probs, 
                         pheno = pheno, 
                         kinship = kinship[[chrom]], 
                         addcovar = covar_matrix,
                         intcovar = interactive_covariate)
  
} else {
  
  # Pull effects
  scan1coef_out <- qtl2::scan1coef(genoprobs = qtl_probs, 
                                   pheno = pheno, 
                                   kinship = kinship[[chrom]], 
                                   addcovar = covar_matrix)
  
  # Pull blups
  scan1blup_out <- qtl2::scan1blup(genoprobs = qtl_probs, 
                                   pheno = pheno, 
                                   kinship = kinship[[chrom]], 
                                   addcovar = covar_matrix)
  
  # Fit effects at peak
  fit1_out <- qtl2::fit1(genoprobs = qtl_peak_probs, 
                         pheno = pheno, 
                         kinship = kinship[[chrom]], 
                         addcovar = covar_matrix)
  
}

# make peak info matrix
qtl_coord <- data.frame(phenotype, pos_start, pos_peak, pos_end)
peaks <- cbind(chrom,qtl_coord,fit1_out$lod,t(fit1_out$coef))
colnames(peaks)[colnames(peaks) == "fit1_out$lod"] <- "LOD"

# Plot effects
png(paste0(phenotype,"_scan1coef.png"))
qtl2::plot_coefCC(x = scan1coef_out, 
                  map = map, 
                  columns = LETTERS[1:8], 
                  scan1_output = scan1out)
dev.off()

# Plot blups
png(paste0(phenotype,"_scan1blup.png"))
qtl2::plot_coefCC(x = scan1blup_out, 
                  map = map, 
                  columns = LETTERS[1:8], 
                  scan1_output = scan1out)
dev.off()

# Save effect files and peak info
save(scan1coef_out, scan1blup_out, fit1_out, 
     file = paste(phenotype, chrom, round(pos_start,2), round(pos_end,2),"qtl_effect_files.RData",sep = "_"))
write.csv(peaks, file = paste(phenotype, chrom, round(pos_start,2), round(pos_end,2),"peaks.csv",sep = "_"), row.names = FALSE, quote = FALSE)


