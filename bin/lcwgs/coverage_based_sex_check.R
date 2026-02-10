#!/usr/bin/env Rscript

################################################################################
# Infer sex of sample from relative coverage on chromosome X.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20250317
################################################################################
library(purrr)
library(dplyr)

# take arguments
args <- commandArgs(trailingOnly = TRUE)

# what chromosome?
metadata <- args[1]
covar <- read.csv(metadata, tryLogical = F)

# sample coverage files
mosdepth_files <- list.files(pattern = "mosdepth.summary.txt")

# calculate relative coverage of all sites on X vs. all sites in the genome
relative_coverage_list <- lapply(mosdepth_files, function(m){
  cov_file = read.table(m, header = T)
  sex_chr = cov_file %>%
    dplyr::filter(chrom %in% c("X","Y","total"))
  x = sex_chr[sex_chr$chrom == "X",]$mean
  y = sex_chr[sex_chr$chrom == "Y",]$mean
  gw = sex_chr[sex_chr$chrom == "total",]$mean
  sex_ratio = x/y
  gw_ratio = x/gw
  id = strsplit(m,"[.]")[[1]][[1]]
  data.frame(id, sex_ratio, gw_ratio)
})
relative_coverage_df <- Reduce(rbind, relative_coverage_list) %>%
  dplyr::left_join(., covar)

# denote sex information
new_covar <- relative_coverage_df %>%
  dplyr::mutate(inferred_sex = dplyr::case_when(gw_ratio > 0.75 ~ "F",
                                                gw_ratio <= 0.75 ~ "M",
                                                is.nan(gw_ratio) | is.na(gw_ratio) ~ sex))

# select variables
sex_checked <- new_covar[,c(colnames(covar),"inferred_sex","gw_ratio","sex_ratio")] %>%
  dplyr::rename(original_sex = sex,
                sex = inferred_sex)
slim_new_covar <- sex_checked[,c(colnames(covar),"original_sex","gw_ratio","sex_ratio")]
  
# write the new covar file
write.csv(slim_new_covar, file = "sex_check_covar.csv", row.names = F, quote = F)


