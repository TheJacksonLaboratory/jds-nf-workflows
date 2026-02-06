#!/usr/bin/env Rscript

################################################################################
# Calculate genotype and allele probabilities from filtered QUILT imputed
# genotype calls.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20260123
################################################################################

library(dplyr)
library(qtl2)

# take arguments
args <- commandArgs(trailingOnly = TRUE)

# what chromosome?
chrom <- args[1]

# sample genotypes
sample_genos <- args[2]

# founder genotypes
founder_genos <- args[3]

# physical map
pmap <- args[4]

# genetic map
gmap <- args[5]

# sample metadata
metadata <- args[6]

# cross type
cross_type <- args[7]

# smoothing function
smooth_function <- args[8]

# smoothing window
smooth_window <- as.numeric(args[9])

# Number of cores from Nextflow
n_cores <- as.numeric(args[10])
message(paste("Using", n_cores, "cores for parallel operations"))

## RUN ##

cat(paste0("Processing chromosome ",chrom,"\n"))
cat("Reading in covariate data...\n")
covar <- read.csv(metadata, tryLogical = F)
if("X" %in% colnames(covar)){
  revised_covar <- covar[,-1]
  write.csv(revised_covar, file = "covar.csv", row.names = F, quote = F)
  covar <- read.csv(metadata, tryLogical = F)
}
if(all(covar$sex == FALSE)){
  covar$sex <- "F"
}
if("original_sex" %in% colnames(covar) & all(covar$original_sex == FALSE)){
  covar$original_sex <- "F"
}
if(cross_type != "do"){
  revised_covar <- covar %>%
    dplyr::select(-original_sex, -gw_ratio, -sex_ratio)
  write.csv(revised_covar, file = "covar.csv", row.names = F, quote = F)
  covar <- read.csv(metadata, tryLogical = F)
}



# Has a cross type been specified?
stopifnot("cross_type" %in% ls(pattern = "cross_type"))

cat("Writing control file for qtl2...\n")
if(cross_type == "genail4" | cross_type == "het3"){
  # write control file for genail4 crosses
  write_control_file(output_file = paste0("chr",chrom,"_control_file.json"),
                     crosstype="genail4",
                     founder_geno_file=founder_genos,
                     founder_geno_transposed=TRUE,
                     gmap_file=gmap,
                     pmap_file=pmap,
                     geno_file=sample_genos,
                     geno_transposed=TRUE,
                     geno_codes=list(A=1, H=2, B=3),
                     sex_covar="sex",
                     sex_codes=list("F"="female", 
                                    "M"="male"),
                     covar_file = metadata,
                     crossinfo_covar = colnames(covar)[!colnames(covar) %in% c("id","sex")],
                     xchr = "X",
                     overwrite = T)
} else if(cross_type == "F1_mut"){
  # write control file for genail4 crosses
  write_control_file(output_file = paste0("chr",chrom,"_control_file.json"),
                     crosstype="genail14",
                     founder_geno_file=founder_genos,
                     founder_geno_transposed=TRUE,
                     gmap_file=gmap,
                     pmap_file=pmap,
                     geno_file=sample_genos,
                     geno_transposed=TRUE,
                     geno_codes=list(A=1, H=2, B=3),
                     sex_covar="sex",
                     sex_codes=list("F"="female", 
                                    "M"="male"),
                     covar_file = metadata,
                     crossinfo_covar = colnames(covar)[!colnames(covar) %in% c("id","sex")],
                     xchr = "X",
                     overwrite = T)
} else if(cross_type == "do"){
  # Write control file for DO crosses
  qtl2::write_control_file(output_file = paste0("chr",chrom,"_control_file.json"),
                           crosstype="do",
                           founder_geno_file=founder_genos,
                           founder_geno_transposed=TRUE,
                           gmap_file=gmap,
                           pmap_file=pmap,
                           geno_file=sample_genos,
                           geno_transposed = TRUE,
                           geno_codes=list(A=1, H=2, B=3),
                           sex_covar = "sex",
                           sex_codes=list("F"="female", 
                                          "M"="male"),
                           covar_file = metadata,
                           crossinfo_covar="gen",
                           xchr = "X",
                           overwrite = T)
} else if(cross_type == "cc"){
  # Write control file for CC
  # write control file for genail4 crosses
  write_control_file(output_file = paste0("chr",chrom,"_control_file.json"),
                     crosstype="genail8",
                     founder_geno_file=founder_genos,
                     founder_geno_transposed=TRUE,
                     gmap_file=gmap,
                     pmap_file=pmap,
                     geno_file=sample_genos,
                     geno_transposed=TRUE,
                     geno_codes=list(A=1, H=2, B=3),
                     sex_covar="sex",
                     sex_codes=list("F"="female", 
                                    "M"="male"),
                     covar_file = metadata,
                     crossinfo_covar = colnames(covar)[!colnames(covar) %in% c("id","sex")],
                     xchr = "X",
                     overwrite = T)
} else if(cross_type == "bxd"){
  qtl2::write_control_file(output_file = paste0("chr",chrom,"_control_file.json"),
                         crosstype="risib",
                         geno_file=sample_genos,
                         geno_transposed=TRUE,
                         geno_codes=list(B=1, D=2),
                         xchr="X",
                         gmap_file=gmap,
                         pmap_file=pmap,
                         crossinfo_file = metadata,
                         crossinfo_codes = c("BxD"=0),
                         alleles=c("B", "D"),
                         overwrite=TRUE)
} else {
  print("Cross type specified does not have a process to make .json file; ending")
}

# test if any sample genotypes from the segment of the chromosome
any_markers <- nrow(read.csv(sample_genos))
if (any_markers == 0) {
  message("No markers in chromosome segment, exiting.")
  quit(save = "no", status = 0, runLast = FALSE)
}

# Load in the cross object
cat("Reading in cross object...\n")
cross <- qtl2::read_cross2(paste0("chr",chrom,"_control_file.json"))

# Drop null markers
cross <- qtl2::drop_nullmarkers(cross)

# subset during testing
# cross <- subset(cross, ind = c(1:3))

# Calculate genotype probs
cat("Calculating genotype probabilities...\n")
pr <- qtl2::calc_genoprob(cross = cross, 
                          map = cross$pmap, 
                          error_prob = 0.002, 
                          cores = n_cores, 
                          quiet = F)

# Smooth genotype probs
cat("Smoothing genotype probabilities...\n")
source(smooth_function)
pr <- smooth_genoprobs(pr, window = smooth_window)

# Save objects
cat("Saving output files...\n")
save(cross, file = paste0("chr_",chrom,"_cross.RData"))
save(pr, file = paste0("chr_",chrom,"_36_state_probs.RData"))
#save(apr, file = paste0("chr_",chrom,"_8_state_probs.RData"))
