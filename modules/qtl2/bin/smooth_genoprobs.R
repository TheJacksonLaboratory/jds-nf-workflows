################################################################################
# Given a qtl2-style genoprobs or allele probs object with markers at a 
# high-density, smooth the probs at a given window size.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2025-04-24
################################################################################

library(accelerometry)
library(qtl2)

# Interpolate the probs object down to a lower marker density, smoothing within
# a given window.
# Arguments:
# probs: qtl2-style probs object (genoprobs or allele probs). Named list in
#        which each element is a 3D array containing the diplotype probs for 
#        one chromosome.
# window: integer that is the width of the smoothing window. Should be an
#         odd number. i.e. 5 means that we will smooth +/- 2 markers around
#         the current marker.
# Returns: qtl2-style probs object containing smoothed probs.
smooth_genoprobs = function(probs, window = 5) {

  # Verify that the objects are the correct class.
  if(!'calc_genoprob' %in% class(probs)) {
  
    stop('probs must be of class calc_genoprob')
  
  } # if(!'calc_genoprob' %in% class(probs))
  
  # Smooth each chromosome.
  new_probs = setNames(as.list(names(probs)), names(probs))

  for(i in seq_along(probs)) {
  
    print(paste('CHR:', names(probs)[i]))
  
    new_probs[[i]] = smooth_one_chr(probs[[i]], window)
  
  } # for(i)
  
  attributes(new_probs) = attributes(probs)
  
  return(new_probs)

} # smooth_genoprobs()


# Smooth one chromosome.
smooth_one_chr = function(pr, win) {

  new_pr = array(0, dim = dim(pr), dimnames = dimnames(pr))

  # Fill in the beginning and the end of the probs.
  half_win  = floor(win / 2)
  start_rng = 1:half_win

  # Looping through the first markers.
  for(i in start_rng) {

    new_pr[,,i] = apply(pr[,,1:i, drop = FALSE], 1:2, mean)

  } # for(i)

  n_markers = dim(pr)[3]
  end_rng = (n_markers - ceiling(half_win)):n_markers

  # Looping through the last markers.
  for(i in end_rng) {

    new_pr[,,i] = apply(pr[,,i:n_markers, drop = FALSE], 1:2, mean)

  } # for(i)

  # Loop through samples.
  rng = (half_win + as.numeric(win %% 2 != 0)):(n_markers - half_win)
  for(i in 1:nrow(pr)) {

      new_pr[i,,rng] = t(apply(pr[i,,], 1, movingaves, window = win))

  } # for(i)

  return(new_pr)

} # smooth_one_chr()


############
# Test code.

#base_dir = '/projects/korstanje-lab/Pureplex/AnalyzedData/TumorStudy_combined/results/quilt/20250421_tumorstudy_combined/2000/geno_probs'
#probs = readRDS(file.path(base_dir, 'complete_alleleprobs.rds'))

#np = smooth_genoprobs(probs, 300)

