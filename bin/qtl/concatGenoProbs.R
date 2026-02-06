#!/usr/bin/env Rscript
library(qtl2)
library(parallel)

################################################################################
# Concatenate project genotype probabilities, convert to allele probabilities, 
# calculate genotyping errors, etc.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20251029
################################################################################ 

# # testing
# test_dir <- "/flashscratch/widmas/HR_QC_outputDir/work/88/c6ee71e4f7c91e94491e64f99873b8"
# setwd(test_dir)

# get genoprobs objects
cross_file_list <- c(t(read.table("crosses.txt")))
cross_list <- lapply(cross_file_list, function(x){
  cross <- readRDS(x)
  return(cross)
})

# combine sample genotypes
new_geno <- list()
for(c in names(cross_list[[1]]$geno)){
  message(paste("Chromosome",c))
  geno <- Reduce(rbind,lapply(cross_list, function(x){
    g <- x$geno[[c]]
    return(g)
  }))
  new_geno[[c]] <- geno
}

# combine covars
new_covar <- Reduce(rbind,lapply(cross_list, function(x){
  cov <- x$covar
  return(cov)
}))

# combine crossinfos
new_crossinfo <- Reduce(rbind,lapply(cross_list, function(x){
  ci <- x$cross_info
  return(ci)
}))

# combine isfemales
new_isfemales <- unlist(lapply(cross_list, function(x) x$is_female))

# skeleton cross to replace with concatenated sample data
cross <- cross_list[[1]]
cross$covar <- new_covar
cross$cross_info <- new_crossinfo
cross$is_female <- new_isfemales
cross$geno <- new_geno

# get genoprobs objects
probs_file_list <- c(t(read.table("probs.txt")))

# read probs in
probs_list <- lapply(probs_file_list, function(x){
  pr <- readRDS(x)
  return(pr)
})

# concatenate
genoprobs <- Reduce(rbind, probs_list)

# clean genotype probabilities
genoprobs <- qtl2::clean_genoprob(genoprobs, cores = parallel::detectCores())

# convert to allele probs
alleleprobs <- qtl2::genoprob_to_alleleprob(probs = genoprobs, 
                                            cores = parallel::detectCores(), 
                                            quiet = F)

# calculate viterbi
m <- qtl2::maxmarg(probs = genoprobs, minprob=0.5, quiet = T)

# calculate kinship matrix
k <- calc_kinship(genoprobs, type = "loco", quiet = FALSE, cores = parallel::detectCores())

# identify missing marker genotypes in the cohort 
percent_missing_marker <- qtl2::n_missing(cross, "marker", "prop")*100
bad_markers <- names(percent_missing_marker[which(percent_missing_marker > 10)])

# calculate genotyping errors
e <- calc_errorlod(cross, genoprobs, cores = parallel::detectCores())

# save
saveRDS(cross$pmap, file = "pmap.rds")
saveRDS(cross$gmap, file = "gmap.rds")
saveRDS(cross, file = "cross.rds")
saveRDS(genoprobs, file = "genoprobs.rds", compress = TRUE)
saveRDS(alleleprobs, file = "alleleprobs.rds", compress = TRUE)
saveRDS(m, file = "maxmarg.rds", compress = TRUE)
saveRDS(e, file = "genotyping_errors.rds", compress = TRUE)
saveRDS(k, file = "kinship.rds", compress = TRUE)
saveRDS(bad_markers, file = "bad_markers.rds")
