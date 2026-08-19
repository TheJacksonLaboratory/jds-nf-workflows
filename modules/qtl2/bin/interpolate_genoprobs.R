################################################################################
# Given two sets of genotype probabilities and thier marker postions,
# interpolate the probabilities of the first set to the second set.
# Daniel Gatti
# dan.gatti@jax.org
# 2022-09-08
################################################################################

library(IRanges)
library(qtl2convert)
library(qtl2)

# Interpolate pr1 onto the positions of mkr2.
# Arguments:
# pr1: 3D numeric array of genoprobs with samples in rows, founders in columns, 
#      & markers in slices.
# mkr1: Named numeric vector of marker position (in bp) in pr1.
# mkr2: Names numeric vector of marker positions (in bp) to interpolate genoprobs to.

interpolate_one_chr = function(pr1, mkr1, mkr2) {
  
  # Verify that the marker names in pr1 and mkr1 are identical and in the same
  # order.
  stopifnot(dimnames(pr1)[[3]] == names(mkr1))
  
  # Verify that the markers are in bp and not Mbp.
  # Chr 1 is 195 Mbp long, so if the minimum value is greater than 200,
  # I'm assuming the marker units are bp. The start of most chromosomes is
  # around 3 Mbp.
  stopifnot(min(mkr1) > 200)
  stopifnot(min(mkr2) > 200)

  # Convert markers to IRanges to used inter-range functions.
  mkr1 = IRanges(start = mkr1, width = 1, names = names(mkr1))
  mkr2 = IRanges(start = mkr2, width = 1, names = names(mkr2))
  
  # The the proximal and distal markers in mkr1 which bracket the markers
  # in mkr2.
  proximal = follow(mkr2,  mkr1)
  distal   = precede(mkr2, mkr1)
  proximal[is.na(proximal)] = 1
  distal[is.na(distal)]     = length(mkr1)
  
  # Ensure indices are within the range
  proximal[proximal > dim(pr1)[3]] = dim(pr1)[3]
  distal[distal     > dim(pr1)[3]] = dim(pr1)[3]
  
  # Get the proximal and distal probs.
  prox_pr = pr1[,,proximal, drop = FALSE]
  dist_pr = pr1[,,distal,   drop = FALSE]
  
  # interpolate values
  mkr_interp = (start(mkr2) - start(mkr1)[proximal]) / (start(mkr1)[distal] - start(mkr1)[proximal])
  mkr_interp[is.infinite(mkr_interp)] = 1
  mkr_interp = array(rep(mkr_interp, each = nrow(pr1) * ncol(pr1)),
                     c(nrow(pr1), ncol(pr1), length(mkr2)),
                     dimnames = list(rownames(pr1), colnames(pr1), names(mkr2)))

  # Create the new probs array.
  tmp = prox_pr + (mkr_interp) * (dist_pr - prox_pr)
  dimnames(tmp)[[3]] = names(mkr2)
  return(tmp)

} # interpolate_one_chr()


# Top-level interpolate function.
# Arguments:
# probs1: qtl2-style genoprobs. Named list containing 3D numerical arrays with
#         samples in rows, founders in columns, and markers in slices. Each
#         list element contains one chromosome and names contain chromosome name. 
# markers1: qtl2-style marker map. Named list containing a named numeric 
#           vector of marker positions for one chromosome. Names of vector 
#           are markers in probs1. Each list element contains one chromosome 
#           and names contain chromosome name.
# markers2: qtl2-style marker map. Named list containing a named numeric 
#           vector of marker positions for one chromosome. Names of vector 
#           are new markers to interpolate to. Each list element contains one 
#           chromosome and names contain chromosome name. probs1 will be 
#           interpolated onto the markers in this map.
interpolate_genoprobs = function(probs1, markers1, markers2) {

  # Perform some minimal error checking. 
  # TBD: More rigorous checks for matching chromosome names in each object.
  if(length(probs1) != length(markers1)) {
    stop('Probs1 and markers1 must be the same length.')
  }

  # If the marker maps are in Mb, convert to bp.
  # Chr 1 is the longest chromsome and it is less than 200 Mb.
  if(max(unlist(markers1) < 200)) {

    markers1 = lapply(markers1, '*', 1e6)

  } # if(max(unlist(markers1) < 200))

  if(max(unlist(markers2) < 200)) {

    markers2 = lapply(markers2, '*', 1e6)

  } # if(max(unlist(markers2) < 200))  
  
  all_chr = names(probs1)

  new_probs1        = vector('list', length(probs1))
  names(new_probs1) = names(probs1)

  for(chr in all_chr) {
      
    new_probs1[[chr]] = interpolate_one_chr(pr1  = probs1[[chr]],
                                            mkr1 = markers1[[chr]],
                                            mkr2 = markers2[[chr]])
    
  } # for(chr)

  # Set the attributes for the probs object.
  attr(new_probs1, 'crosstype')    = attr(probs1, 'crosstype')
  attr(new_probs1, 'is_x_chr')     = attr(probs1, 'is_x_chr')
  attr(new_probs1, 'alleles')      = attr(probs1, 'alleles')
  attr(new_probs1, 'alleleprobs')  = attr(probs1, 'alleleprobs')
  class(new_probs1)                = class(probs1)
  
  return(new_probs1)
  
} # interpolate_genoprobs()


