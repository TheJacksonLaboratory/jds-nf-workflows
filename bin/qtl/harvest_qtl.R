#!/usr/bin/env Rscript

################################################################################
# Gather and summarize QTL results
# - Get thresholds, read in LODs, make new plots with thresholds on them.
# - Harvest QTL peaks, get 95% Bayesian credible intervals, and write out table
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20250626
################################################################################

library(qtl2)

# Test files
scan1_files <- list.files(pattern = "scan1out.rds")
perm_files <- list.files(pattern = "scan1perms.rds")
cross_files <- list.files(pattern = "cross.rds")

# Grab phenotype names
phenotypes <- unlist(lapply(perm_files, function(x){
  pheno <- gsub(pattern = "_scan1perms.rds",replacement = "", basename(x))
  return(pheno)
}))

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
  peaks[[i]] <-  qtl2::find_peaks(scan1_output = scan1, 
                                  map = cross$pmap,
                                  threshold = summary(perms)[[1]], drop = 3, peakdrop = 3,
                                  thresholdX = summary(perms)[[2]], dropX = 3, peakdropX = 3,
                                  expand2markers = TRUE,
                                  sort_by = "pos")
  
}
QTL_table <- Reduce(rbind, peaks)
write.csv(QTL_table, file = "peaks.csv", row.names = F, quote = F)

