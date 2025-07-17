#!/usr/bin/env Rscript

################################################################################
# Perform data quality control for QTL analysis
# - Read in all of the data
# - Verify that phenotypes, covar, and probs have the same samples
# - Subset to the common samples and order the samples the same way in all objects.
# - Verify that map and probs have the same markers 
# - Verify that the covar info/transformation/mapping model matches the phenotypes and covariates that we have.
# - Transform the phenotypes, if requested.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20250716
################################################################################

library(qtl2)
library(dplyr)

# Read in the data
args <- commandArgs(trailingOnly = TRUE)

# # Inputs (for testing):
# baseDir           <- "/flashscratch/widmas/qtl_mapping_outputDir/work/3d/52883ff42ddafa870f93cadf526e8d/"
# setwd(baseDir)
# pheno_file        <- "attie_test_pheno.csv"
# covar_file        <- "attie_500_test_covar.csv"
# genoprobs_file    <- "attie_500_test_genoprobs.rds"
# alleleprobs_file  <- "attie_500_test_alleleprobs.rds"
# cross_file        <- "attie_500_test_cross.rds"
# kinship_file      <- "attie_500_test_kinship.rds"
# pheno_file        <- "attie_test_pheno.csv"
# covar_info_file    <- "attie_test_transform.csv"
# outdir            <- "/flashscratch/widmas/qtl_mapping_qc"
# dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
covar_file              <- args[1]
cross_file              <- args[2]
genoprobs_file          <- args[3]
alleleprobs_file        <- args[4]
kinship_file            <- args[5]
pheno_file              <- args[6]
covar_info_file         <- args[7]

# Read in the phenotype data
pheno <- read.csv(pheno_file)
covar_info <- read.csv(covar_info_file)
transform <- covar_info %>%
  dplyr::distinct(phenotype, transformation)

# Fix phenotype columns
colnames(pheno) <- gsub(" ","_",colnames(pheno))
transform$phenotype <- gsub(" ","_",transform$phenotype)
stopifnot(all(transform$phenotype %in% colnames(pheno)))

# Read in the covariate data
covar <- read.csv(covar_file)

common_samples <- intersect(pheno$id, covar$id)
# Test if all the samples in phenotype file are in covar file
if (!all(pheno$id %in% covar$id)) {
  print("Not all samples in phenotype file are in covariate file.")
  pheno <- pheno[pheno$id %in% common_samples, ]
}
# Test if all the samples in covar file are in phenotype file
if (!all(covar$id %in% pheno$id)) {
  print("Not all samples in covariate file are in phenotype file.")
  covar <- covar[covar$id %in% common_samples, ]
}

# Order the samples the same way in both objects
pheno <- pheno[order(pheno$id), ]
covar <- covar[order(covar$id), ]
stopifnot(pheno$id == covar$id)

# Transform the phenotypes according to the transform file
transformed_phenos <- list()
for(i in 1:nrow(transform)){
  t <- transform[i,]$transformation
  p <- transform[i,]$phenotype
  
  if(is.na(t)){
    
    untrans_pheno <- data.frame(pheno[,which(colnames(pheno) == p)])
    colnames(untrans_pheno) <- colnames(pheno)[which(colnames(pheno) == p)]
    transformed_phenos[[i]] <- untrans_pheno
    
  } else if(t == "log"){
    
    trans_pheno <- data.frame(log10(pheno[,which(colnames(pheno) == p)]))
    colnames(trans_pheno) <- paste0("log_",p)
    transformed_phenos[[i]] <- trans_pheno
    covar_info$phenotype[which(covar_info$phenotype == p & covar_info$transformation == t)] <- paste0("log_",p)
    
  } else if(t == "log1p"){
    
    trans_pheno <- data.frame(log1p(pheno[,which(colnames(pheno) == p)]))
    colnames(trans_pheno) <- paste0("log1p_",p)
    transformed_phenos[[i]] <- trans_pheno
    covar_info$phenotype[which(covar_info$phenotype == p & covar_info$transformation == t)] <- paste0("log1p_",p)
    
  } else if(t == "sqrt"){
    
    trans_pheno <- data.frame(sqrt(pheno[,which(colnames(pheno) == p)]))
    colnames(trans_pheno) <- paste0("sqrt_",p)
    transformed_phenos[[i]] <- trans_pheno
    covar_info$phenotype[which(covar_info$phenotype == p & covar_info$transformation == t)] <- paste0("sqrt_",p)
    
  } else if(t == "rankZ"){
    
  } else if(t == "exp"){
    
    trans_pheno <- data.frame(exp(pheno[,which(colnames(pheno) == p)]))
    colnames(trans_pheno) <- paste0("e_",p)
    transformed_phenos[[i]] <- trans_pheno
    covar_info$phenotype[which(covar_info$phenotype == p & covar_info$transformation == t)] <- paste0("e_",p)
    
  } else {
    
    stop("One of log, log1p, sqrt, rankZ, exp tranformations, or NA (no transformation) not specified; quitting")
  }
  
}
transformed_phenos <- Reduce(cbind, transformed_phenos)
pheno <- cbind(pheno, transformed_phenos) %>%
  dplyr::select(id, colnames(transformed_phenos))


# Read in the genotype probabilities
genoprobs <- readRDS(genoprobs_file)
genoprobs <- subset(genoprobs, ind = covar$id)

# Read in the allele probabilities
alleleprobs <- readRDS(alleleprobs_file)
alleleprobs <- subset(alleleprobs, ind = covar$id)

# Read in the cross object
cross <- readRDS(cross_file)
cross <- subset(cross, ind = covar$id)

# Read in the kinship matrix
kinship <- readRDS(kinship_file)
stopifnot(unique(unlist(lapply(kinship, dim))) == nrow(covar))
kinship <- lapply(kinship, function(x){
  new_k <- x[covar$id,covar$id]
  str(x)
  str(new_k)
  attr(new_k,"n_pos") <- attr(x,"n_pos")
  return(new_k)
})

# Verify samples are in the same order in the probs objects
stopifnot(dimnames(genoprobs[[1]])[[1]] == dimnames(alleleprobs[[1]])[[1]])

# Verify markers are the same in each of the probs objects
for(i in names(genoprobs)){
  stopifnot(dimnames(genoprobs[[i]])[[1]] == dimnames(alleleprobs[[i]])[[1]])
  stopifnot(dimnames(genoprobs[[i]])[[3]] == dimnames(alleleprobs[[i]])[[3]])
}

# make phenotype files for each phenotype
pheno_cols <- colnames(pheno)[!colnames(pheno) %in% c("id","sex","gen")]
for(i in pheno_cols){
  write.csv(pheno[,c("id",i)], file = paste0(i,"_pheno.csv"), row.names = FALSE, quote = FALSE)
}

for(i in pheno_cols){
  filtered_covar_info <- covar_info %>%
    dplyr::filter(phenotype == i) %>%
    dplyr::select(-transformation)
  write.csv(filtered_covar_info, file = paste0(i,"_covar_info.csv"), row.names = FALSE, quote = FALSE)
}

# Save the updated files
saveRDS(genoprobs, file = "pr.rds")
saveRDS(alleleprobs, file = "apr.rds")
saveRDS(cross, file = "cross.rds")
saveRDS(kinship, file = "kinship.rds")
write.csv(covar, file = "covar.csv", row.names = FALSE, quote = FALSE)
