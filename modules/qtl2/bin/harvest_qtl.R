#!/usr/bin/env Rscript

################################################################################
# Gather and summarize QTL results
# - Get thresholds, read in LODs, make new plots with thresholds on them.
# - Harvest QTL peaks, get 95% Bayesian credible intervals, and write out table
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20260311
################################################################################

library(qtl2)

# Inputs
scan1_files <- list.files(pattern = "scan1out.rds")
perm_files <- list.files(pattern = "scan1perms.txt")
map_files <- list.files(pattern = "map.rds")

# Grab phenotype names
phenotypes <- unlist(lapply(scan1_files, function(x){
  pheno <- gsub(pattern = "_scan1out.rds",replacement = "", basename(x))
  return(pheno)
}))

# Pull significant QTLs and Find peaks
peaks <- list()
for(i in phenotypes){
  
  # read in data
  
  scan1 <- readRDS(paste0(i,"_scan1out.rds"))
  perms <- read.table(paste0(i,"_scan1perms.txt"), header = T)
  map   <- readRDS( paste0(i,"_map.rds"))
  
  # plot
  png(paste0(i,"_scan1_thresh.png"))
  qtl2::plot_scan1(scan1, map, main = i)
  qtl2::add_threshold(map, 
                      thresholdA = perms$perm, 
                      thresholdX = perms$perm, col = "red")
  dev.off()
  
  # find peaks
  peaks[[i]] <-  qtl2::find_peaks(scan1_output = scan1, 
                                  map = map,
                                  threshold = perms$perm, drop = 3, peakdrop = 3,
                                  thresholdX = perms$perm, dropX = 3, peakdropX = 3,
                                  expand2markers = TRUE,
                                  sort_by = "pos")
  
}
QTL_table <- Reduce(rbind, peaks)
write.csv(QTL_table, file = "peaks.csv", row.names = F, quote = F)

