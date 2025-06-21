#!/usr/bin/env Rscript

################################################################################
# Gather and summarize QTL results
# - Get thresholds, read in LODs, make new plots with thresholds on them.
# - Harvest QTL peaks, get 95% Bayesian credible intervals, and write out table
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20250619
################################################################################

library(qtl2)

# setwd("/flashscratch/widmas/outputDir/work/78/9df6e92e90bfb09251cb28f303da77")

# Test files
scan1_files <- list.files(pattern = "scan1out.rds")
perm_files <- list.files(pattern = "scan1perms.rds")
cross_files <- list.files(pattern = "cross.rds")

# Grab phenotype names
phenotypes <- unlist(lapply(perm_files, function(x){
  pheno <- gsub(pattern = "_scan1perms.rds",replacement = "", basename(x))
  return(pheno)
}))

# Plot QTL scans with permutation thresholds
for(i in phenotypes){
  scan1 <- readRDS(scan1_files[grep(i, scan1_files)])
  perms <- readRDS(perm_files[grep(i, perm_files)])
  cross <- readRDS(cross_files[grep(i, cross_files)])
  thresh <- as.numeric()
  png(paste0(i,"_scan1_thresh.png"))
  qtl2::plot_scan1(scan1, cross$pmap, main = i)
  qtl2::add_threshold(cross$pmap, 
                      thresholdA = summary(perms)[[1]], 
                      thresholdX = summary(perms)[[2]], col = "red")
  dev.off()
}

# Pull significant QTLs and Find peaks
peaks <- list()
for(i in phenotypes){
  
  # read in data
  scan1 <- readRDS(scan1_files[grep(i, scan1_files)])
  perms <- readRDS(perm_files[grep(i, perm_files)])
  cross <- readRDS(cross_files[grep(i, cross_files)])
  
  # plot
  png(paste0(i,"_scan1_thresh.png"))
  qtl2::plot_scan1(scan1, cross$pmap, main = i)
  qtl2::add_threshold(cross$pmap, 
                      thresholdA = summary(perms)[[1]], 
                      thresholdX = summary(perms)[[2]], col = "red")
  dev.off()
  
  # find peaks
  peaks[[i]] <- find_peaks(scan1_output = scan1, 
                           map = cross$pmap, 
                           threshold  = summary(perms)[[1]], prob = 0.95, 
                           thresholdX = summary(perms)[[2]], probX = 0.95,
                           sort_by = "pos",
                           cores = 0)
  
}
QTL_table <- Reduce(rbind, peaks)
write.csv(QTL_table, file = "peaks.csv")

