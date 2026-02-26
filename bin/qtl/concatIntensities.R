#!/usr/bin/env Rscript
library(fst)
library(parallel)
library(dplyr)
library(purrr)

################################################################################
# Concatenate project intensity files.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20250521
################################################################################ 

# testing
# test_dir <- "/flashscratch/widmas/HR_QC_outputDir/work/e0/fcecb9bbc5d0c894e96155b0caaa00"
# setwd(test_dir)

# read in intensities and make longer to easily summarize duplicate stats
parseSexChrIntsLong <- function(int_file){
  # extract the finalreport file hash
  hash <- strsplit(x = int_file, split = "_")[[1]][[2]]
  
  # make the intensities long
  long_int <- read.csv(int_file, skip = 3, check.names = F) %>%
    dplyr::mutate(hash = hash) %>%
    tidyr::pivot_longer(-c(marker, hash), names_to = "sample", values_to = "int")
  
  return(long_int)
  
}


# get chr X intensities objects
x_files <- list.files(pattern = "chrX")
x_int_list_long <- purrr::map(x_files, parseSexChrIntsLong)
x_int_df_long <- Reduce(rbind, x_int_list_long)

# get samples that were genotyped multiple times
sample_counts <- x_int_df_long %>%
  dplyr::distinct(sample, hash) %>%
  dplyr::group_by(sample) %>%
  dplyr::count()

# find the run with the cleanest intensities
dedup_samples <- sample_counts %>%
  dplyr::left_join(., x_int_df_long) %>%
  dplyr::select(-n) %>%
  dplyr::arrange(marker) %>%
  dplyr::group_by(sample, hash) %>% ## count the NA intensities and estimate intensity variance
  dplyr::summarize(int_variance = var(int, na.rm = TRUE),
                   na_count = sum(is.na(int)),
                   .groups = 'drop') %>%
  dplyr::group_by(sample) %>%      ## get the sample with the fewest NAs
  dplyr::slice_min(order_by = na_count, n = 1) %>% ## if there are ties get the sample with the lowest intensity variance
  dplyr::slice_min(order_by = int_variance, n = 1) %>%
  dplyr::select(sample, hash)

# join back to the intensities for the hash of interest
x_int_df <- dedup_samples %>%
  dplyr::left_join(., x_int_df_long) %>%
  dplyr::select(-hash) %>%
  tidyr::pivot_wider(names_from = sample, values_from = int) %>%
  data.frame(., check.names = F)

# get chr Y intensities objects
y_files <- list.files(pattern = "chrY")
y_int_list_long <- purrr::map(y_files, parseSexChrIntsLong)
y_int_df_long <- Reduce(rbind, y_int_list_long)
y_int_df <- dedup_samples %>%
  dplyr::left_join(., y_int_df_long) %>%
  dplyr::select(-hash) %>%
  tidyr::pivot_wider(names_from = sample, values_from = int) %>%
  data.frame(., check.names = F)

# .fst intensity files
intensity_files <- list.files(pattern = "fst")
filterFsts <- function(int_file){
  # extract the finalreport file hash
  hash <- gsub(pattern = ".fst", replacement = "",strsplit(x = int_file, split = "_")[[1]][[2]])
  
  # read all the intensities
  int <- fst::read.fst(int_file)
  hash_clean_sample_ind <- which(colnames(int) %in% dedup_samples$sample[which(dedup_samples$hash == hash)])
  filtered_int <- int[,c(1:2,hash_clean_sample_ind)]
  
  return(filtered_int)
}
intensity_list <- purrr::map(intensity_files, filterFsts)
int_df <- Reduce(left_join, intensity_list)

# save
write.csv(x_int_df, file = "chrX_intensities.csv", row.names = F, quote = F)
write.csv(y_int_df, file = "chrY_intensities.csv", row.names = F, quote = F)
fst::write.fst(int_df, path = "all_intensities.fst")
write.csv(dedup_samples, file = "dedup_samples.csv", row.names = F, quote = F)
