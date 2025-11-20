#!/usr/bin/env Rscript

################################################################################
# Concatenate genotype probabilities and physical maps from QUILT.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20250327
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
# chroms <- c(as.character(c(18:19)),"X")

# genotype prob objects
genoprobs <- paste0("chr_",chroms,"_36_state_probs.RData")
# genoprobs <- unlist(lapply(genoprobs, function(x) list.files(rundir, pattern = x, recursive = T, full.names = T)))
print(genoprobs)

# read probs in
probs <- vector(mode = "list", length = length(genoprobs))
names(probs) <- chroms
for(i in 1:length(names(probs))){
  message(genoprobs[i])
  load(genoprobs[i])
  probs[[names(probs)[i]]] <- pr[[1]]
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
# pmaps <- unlist(lapply(pmaps, function(x) list.files(rundir, pattern = x, recursive = T, full.names = T)))
new_pmaps <- vector(mode = "list", length = length(pmaps))
names(new_pmaps) <- chroms
for(i in 1:length(names(new_pmaps))){
  load(pmaps[i])
  new_pmaps[[names(new_pmaps)[i]]] <- cross$pmap[[1]]
}
message("Saving genotype probabilities")
saveRDS(object = probs, file = "complete_genoprobs.rds")
saveRDS(object = new_pmaps, file = "complete_pmap.rds")


# make allele probs object
message("Generating allele probabilities")
pr <- probs
rm(probs)
apr <- qtl2::genoprob_to_alleleprob(probs = pr, quiet = F, cores = parallel::detectCores()/1.2)
message("Saving allele probabilities")
saveRDS(object = apr, file = "complete_alleleprobs.rds")

# interpolate
interp_script <- args[3]
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
pr_interp <- interpolate_genoprobs(probs1 = pr,
                                    markers1 = interp_pmap,
                                    markers2 = grid_map)
apr_interp <- interpolate_genoprobs(probs1 = apr,
                                    markers1 = interp_pmap,
                                    markers2 = grid_map)
# divide to get it back to Mb
grid_map <- lapply(grid_map, function(x) x/1e6)

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