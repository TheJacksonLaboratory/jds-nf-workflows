#!/usr/bin/env Rscript

################################################################################
# Run QUILT on bam files
#
# Sam Widmayer
# samuel.widmayer@jax.org
################################################################################


# Load required packages
library(QUILT)

# Set arguments
args <- commandArgs(trailingOnly = TRUE)

#args[1] = bamlist (from channel)
#args[2] = chromosome number (from channel)
#args[3] = covar file
#args[4] = cross type
#args[5] = reference haplotype file
#args[6] = reference sample file
#args[7] = reference legend file
#args[8] = bin shuffle radius

# Input files
# Path to original bamlist
mouse_bamlist <- args[1]

# Chromosome to be analyzed
mouse_chr <- paste0(args[2])

# Path to qtl2 covar file
meta_file <- args[3]

# Read covar file
covar <- read.csv(meta_file)

# Cross type
cross_type <- args[4]

# Reference files
hap   <- args[5]
samp  <- args[6]
leg   <- args[7]
rad   <- args[8]

# Region data
options(scipen = 99999999)
start <- as.integer(args[9])
end   <- as.integer(args[10])
print(args)
# str(args[9])
# str(args[10])
# str(start)
# str(end)

# Take cross type and determine how QUILT should be executed
if(cross_type == "do" | cross_type == "cc" | cross_type == "het3" | cross_type == "F1_mut"){
  
  # DO, CC, and het3 covar files require generation estimates, so can use this for QUILT
  
  # Use median generation from covar file as the nGen parameter in QUILT
  # kill if no generation is included in covar file
  stopifnot("gen" %in% colnames(covar))
  
  # kill if there generation covariate isn't numeric or integer
  stopifnot(is.numeric(covar$gen) || is.integer(covar$gen))
  stopifnot(length(covar$gen[!is.na(covar$gen)]) > 0)
  covar_nGen <- median(covar$gen[!is.na(covar$gen)])
  
  QUILT::QUILT(chr = mouse_chr,
               regionStart = start,
               regionEnd = end,
               buffer = 10000,
               bamlist = mouse_bamlist,
               outputdir = paste0(getwd(), "/"),
               reference_haplotype_file = hap,
               reference_sample_file = samp,
               reference_legend_file = leg,
               shuffle_bin_radius = as.numeric(rad),
               nGen = covar_nGen,
               addOptimalHapsToVCF=TRUE,
               record_interim_dosages=TRUE,
               save_prepared_reference=TRUE)
  
} else if(cross_type == "bxd"){
  
  # BXD covar files do *not* require generation estimates
  # using roughly 20 generations to fully inbred for now
  covar_nGen <- 20
  
  QUILT::QUILT(chr = mouse_chr,
               regionStart = start, # remove when not testing
               regionEnd = end, # remove when not testing
               buffer = 10000, # remove when not testing
               bamlist = mouse_bamlist,
               outputdir = paste0(getwd(), "/"),
               reference_haplotype_file = hap,
               reference_sample_file = samp,
               reference_legend_file = leg,
               shuffle_bin_radius = as.numeric(rad),
               nGen = covar_nGen,
               addOptimalHapsToVCF=TRUE,
               record_interim_dosages=TRUE,
               save_prepared_reference=TRUE)
  
} else {
  
  print("Cross type specified has no QUILT implementation at this time.")
  
}


