#!/usr/bin/env Rscript

################################################################################
# Concatenate genotype probabilities and physical maps from QUILT.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20260120
################################################################################

library(qtl2)
library(tidyr)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)

# cross type
cross_type <- args[1]
# cross_type = "do"

# interpolation gridfile
interp_gridfile <- args[2]

# interpolation script
interp_script <- args[3]

# Number of cores from Nextflow
n_cores <- as.numeric(args[4])
message(paste("Using", n_cores, "cores for parallel operations"))

# cross type
# cross_type <- args[1]
if(cross_type == "het3" | cross_type == "cc" | cross_type == "genail4"){
  cross_type <- "genail4"
} else if(cross_type == "bxd"){
  cross_type <- "risib"
} else if(cross_type == "F1_mut"){
  cross_type <- "genail14"
}else if(cross_type == "do"){
  cross_type <- "do"
} else {
  "No clue of cross type!"
}

# chromosomes
chroms <- c(as.character(c(1:19)),"X")

# genotype prob objects
genoprobs <- paste0("chr_",chroms,"_36_state_probs.RData")
print(genoprobs)

# read probs in
probs <- vector(mode = "list", length = length(genoprobs))
names(probs) <- chroms
for(i in 1:length(names(probs))){
  message(genoprobs[i])
  load(genoprobs[i])
  probs[[names(probs)[i]]] <- pr[[1]]
  rm(pr)  # Remove intermediate object immediately
  gc()    # Force garbage collection
}

# assign attributes
message("Assigning genoprobs attributes...")
attr(probs, "crosstype") <- cross_type
attr(probs, "is_x_chr") <- c(rep(FALSE,19),TRUE)
attr(probs, "alleleprobs") <- FALSE
class(probs) <- c("calc_genoprob", "list")
attr(probs, "alleles") <- unique(unlist(lapply(dimnames(probs[[1]])[[2]],
                                               function(x) strsplit(x, split = "")[[1]])))

# combine physical maps
message("Combining physical maps...")
pmaps <-  paste0("chr_",chroms,"_cross.RData")
new_pmaps <- vector(mode = "list", length = length(pmaps))
names(new_pmaps) <- chroms
for(i in 1:length(names(new_pmaps))){
  load(pmaps[i])
  new_pmaps[[names(new_pmaps)[i]]] <- cross$pmap[[1]]
  rm(cross)  # Remove intermediate object immediately
  gc()       # Force garbage collection
}
is_x_chr <- vector("logical",length = length(new_pmaps))
names(is_x_chr) <- names(new_pmaps)
is_x_chr[which(names(is_x_chr) == "X")] <- TRUE
attr(new_pmaps, "is_x_chr") <- is_x_chr

message("Saving genotype probabilities")
saveRDS(object = probs, file = "complete_genoprobs.rds")
saveRDS(object = new_pmaps, file = "complete_pmap.rds")


# make allele probs object
message("Generating allele probabilities")
apr <- qtl2::genoprob_to_alleleprob(probs = probs, quiet = F, cores = n_cores)
message("Saving allele probabilities")
saveRDS(object = apr, file = "complete_alleleprobs.rds")

# interpolate
source(interp_script)
grid <- read.csv(interp_gridfile)

# filter things
filtered_grid <- grid[grid$chr %in% chroms,]
grid_map <- lapply(unique(filtered_grid$chr), function(x){
  m <- filtered_grid[filtered_grid$chr == x,]$pos*1e6
  names(m) <- filtered_grid[filtered_grid$chr == x,]$marker
  return(m)
})
names(grid_map) <- chroms

# multiply the map to make the ranges work
interp_pmap <- lapply(new_pmaps, function(x) x*1e6)

message("Interpolating genotype probabilities...")
pr_interp <- interpolate_genoprobs(probs1 = probs,
                                    markers1 = interp_pmap,
                                    markers2 = grid_map)
rm(probs)  # Remove original probs after interpolation
gc()

message("Interpolating allele probabilities...")
apr_interp <- interpolate_genoprobs(probs1 = apr,
                                    markers1 = interp_pmap,
                                    markers2 = grid_map)
rm(apr, interp_pmap)  # Clean up
gc()

# divide to get it back to Mb
grid_map <- lapply(grid_map, function(x) x/1e6)
is_x_chr <- vector("logical",length = length(grid_map))
names(is_x_chr) <- names(grid_map)
is_x_chr[which(names(is_x_chr) == "X")] <- TRUE
attr(grid_map, "is_x_chr") <- is_x_chr

# assign interpolated attributes
message("Assigning interpolated genoprob attributes...")
attr(pr_interp, "crosstype") <- cross_type
attr(pr_interp, "is_x_chr") <- c(rep(FALSE,19),TRUE)
attr(pr_interp, "alleleprobs") <- FALSE
class(pr_interp) <- c("calc_genoprob", "list")
attr(pr_interp, "alleles") <- unique(unlist(lapply(dimnames(pr_interp[[1]])[[2]],
                                               function(x) strsplit(x, split = "")[[1]])))

# assign interpolated attributes
message("Assigning interpolated alleleprob attributes...")
attr(apr_interp, "crosstype") <- cross_type
attr(apr_interp, "is_x_chr") <- c(rep(FALSE,19),TRUE)
attr(apr_interp, "alleleprobs") <- TRUE
class(apr_interp) <- c("calc_genoprob", "list")
attr(apr_interp, "alleles") <- unique(unlist(lapply(dimnames(apr_interp[[1]])[[2]],
                                               function(x) strsplit(x, split = "")[[1]])))

# save everything
message("Saving interpolated allele probabilities and pmap")
saveRDS(object = pr_interp, file = "interp_250k_genoprobs.rds")
saveRDS(object = apr_interp, file = "interp_250k_alleleprobs.rds")
saveRDS(object = grid_map, file = "grid_pmap.rds")

# Make kinship matrix
message("Calculating kinship matrix...")
kinship <- qtl2::calc_kinship(probs = apr_interp, 
                              type = "loco", 
                              use_allele_probs = TRUE, 
                              cores = n_cores)
message("Saving kinship matrix...")
saveRDS(object = kinship, file = "interp_250k_kinship_loco.rds")
rm(kinship)  # Clean up after saving
gc()

# Make viterbi object
message("Calculating maxmarg...")
m <- qtl2::maxmarg(probs = pr_interp,
                   minprob = 0.5, 
                   quiet = TRUE, 
                   cores = n_cores)
message("Saving viterbi...")
saveRDS(object = m, file = "interp_250k_viterbi.rds")