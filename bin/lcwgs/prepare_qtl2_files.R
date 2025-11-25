#!/usr/bin/env Rscript

################################################################################
# Gather genotypes, markers, & covariates and format for qtl2::read_cross().
# GRCm39.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20251124
################################################################################

library(qtl2convert)
library(qtl2)
library(stringr)
library(tidyr)
library(dplyr)


# ##### VARIABLES #####

# Arguments:
# founder_file: path to founder VCF.
# sample_file:  path to sample VCF produced by QUILT.
# meta_file:    path to sample metadata file.
# marker_file:  path to QUILT marker file with bp and cM values.
# cross_type:   cross type according to wrt qtl2 preferences

args = commandArgs(trailingOnly = TRUE)

# Founder genotypes and marker positions
founder_file = args[1]

# Sample genotypes from QUILT.
sample_file = args[2]

# Sample metadata file.
meta_file = args[3]

# Cross type
cross_type = args[4]

# Marker map.
marker_file = args[5]

# Chromosome
chr = args[6]

# Grid file
gridfile <- args[7]

# ##### MAIN #####

message("Read in metadata")
meta <- read.csv(meta_file)
if(all(meta$sex == FALSE)){
  meta$sex <- "F"
}

if("original_sex" %in% colnames(meta) & all(meta$original_sex == FALSE)){
  meta$original_sex <- "F"
}

# can't have duplicate sample ids in the metadata for qtl2
meta <- meta[which(!duplicated(meta$id)),]

# can't have any missing values in the covariate file and/or cross info columns
meta <- meta[complete.cases(meta),]

# Read in founder genotypes.
message("Read in founder genotypes")
founder_gt = read.delim(founder_file, check.names = F)
fg_mkrs = paste(founder_gt$CHROM,founder_gt$POS,founder_gt$REF,founder_gt$ALT, sep = "_")
rownames(founder_gt) <- fg_mkrs
quilt_variants <- nrow(founder_gt)

# Read sample genotypes
message("Reading in sample genotypes")
sample_gt <- read.delim(sample_file, check.names = F)
sample_gt <- sample_gt[order(sample_gt$POS),]
sg_mkrs = paste(sample_gt$CHROM,sample_gt$POS,sample_gt$REF,sample_gt$ALT, sep = "_")
rownames(sample_gt) <- sg_mkrs
stopifnot(rownames(sample_gt) == rownames(founder_gt))
stopifnot("id" %in% colnames(meta))
covar_sample_inds <- c()

# get metadata
sample_gt_meta <- sample_gt %>%
  dplyr::select(CHROM, POS, REF, ALT, HWE, INFO_SCORE)

# get genos to fix names
sample_gt_genos <- sample_gt %>%
  dplyr::select(-c(CHROM, POS, REF, ALT, HWE, INFO_SCORE))
stopifnot(rownames(sample_gt_meta) == rownames(sample_gt_genos))

# attempt a clean join; this is only possible when exact sample names are known
if(!all(meta$id %in% colnames(sample_gt_genos)) && !all(colnames(sample_gt_genos) %in% meta$id)){

  message("At least one sample name as supplied in metadata doesn't match sequencing data.")
  message("Searching sample genotypes for metadata sample names...")

  # link the sample IDs
  # assuming that the covar file has some element of the library name in it
  covar_sample_ids <- lapply(colnames(sample_gt_genos), function(x){

    # try a symbol split
    first_split <- stringr::str_split(string = x, pattern = "[:punct:]")[[1]]
    sample_search_1 <- lapply(first_split, function(y){
      putative_id <- meta$id[meta$id %in% y]
      if(length(putative_id) > 0){
        return(putative_id)
      }
    })
    result <- unlist(sample_search_1[which(lapply(sample_search_1, function(x) is.null(x)) == FALSE)])

    # try a dot split
    dot_split <- stringr::str_split(string = x, pattern = "[.]")[[1]]
    sample_search_2 <- lapply(dot_split, function(y){
      putative_id <- meta$id[meta$id %in% y]
      if(length(putative_id) > 0){
        return(putative_id)
      }
    })
    result2 <- unlist(sample_search_2[which(lapply(sample_search_2, function(x) is.null(x)) == FALSE)])

    # try a simple grep
    sample_search_3 <- which(lapply(meta$id, function(y){grep(pattern = paste0("*",y), x = x)}) == 1)
    result3 <- meta$id[sample_search_3]

    if(!is.null(result) && length(result) == 1){
      #print(result)
      return(result)
    } else if (!is.null(result2) && length(result2) == 1){
      #print(result2)
      return(result2)
    } else if (!is.null(result3) && length(result3) == 1){
      #print(result3)
      return(result3)
    } else {
      return(NULL)
    }
  })

  # get rid of samples from quilt vcf not present in the metadata
  # remove the null sample names from the covar file ids
  quilt_samples_absent_from_meta <- unlist(lapply(covar_sample_ids, is.null))
  if(any(quilt_samples_absent_from_meta)){
    covar_sample_ids <- covar_sample_ids[-which(quilt_samples_absent_from_meta)]
    sample_gt_genos <- sample_gt_genos[,-which(quilt_samples_absent_from_meta)]
  }
  # assign the sample names from the metadata present in both files to the sample genotypes
  colnames(sample_gt_genos) <- unlist(covar_sample_ids)

} else {
  message("All sample names as supplied in metadata match sequencing data.")
  sample_gt_genos <- sample_gt_genos
}

# bind things back together
sample_gt_renamed <- cbind(sample_gt_meta, sample_gt_genos)

# recode sample genotypes
message("Recoding sample genotypes.")
recoded_sample_genos <- apply(sample_gt_renamed, 1, function(x){
  REF = x[3]
  ALT = x[4]
  g = as.numeric(gsub(pattern = "[|]",replacement = "",x[-c(1:6)]))
  g[g == 0] <- paste0(REF,REF)
  g[g == 11] <- paste0(ALT,ALT)
  g[g %in% c(10,1)] <- paste0(REF, ALT)
  return(g)
})
recoded_sample_genos <- t(recoded_sample_genos)
colnames(recoded_sample_genos) <- colnames(sample_gt_genos)
sample_gt_renamed <- cbind(sample_gt_meta, recoded_sample_genos)

# recode founder genotypes
message("Recoding founder genotypes.")
recoded_founder_genos <- apply(founder_gt, 1, function(x){
  REF = x[3]
  ALT = x[4]
  g = as.numeric(gsub(pattern = "[|]",replacement = "",x[-c(1:4)]))
  g[g == 0] <- paste0(REF,REF)
  g[g == 11] <- paste0(ALT,ALT)
  g[g %in% c(10,1)] <- paste0(REF, ALT)
  return(g)
})
recoded_founder_genos <- t(recoded_founder_genos)
colnames(recoded_founder_genos) <- colnames(founder_gt)[-c(1:4)]
founder_gt_renamed <- cbind(founder_gt[,1:4], recoded_founder_genos)

# filter metadata down to samples in the genotype file
meta <- meta[which(meta$id %in% colnames(sample_gt_renamed)),]

# at this point, metadata sample size should match the quilt sample size
stopifnot(ncol(sample_gt_genos) == nrow(meta))

# how many sites deviate from HWE?
if(chr != "X"){
  if(cross_type == "do"){

    # filtering by HWE
    sample_gt_renamed <- sample_gt_renamed[which(sample_gt_renamed$HWE > 0.05),]

  } else if(cross_type == "bxd"){
    print("BXD strains; sites not expected to adhere to HWE")
    print("Skipping HWE filter")
    # paste0(round((table(sample_gt_renamed$HWE < 0.05)[[2]]/quilt_variants*100),2),"%")
  } else if(cross_type == "cc"){
    print("CC strains; sites not expected to adhere to HWE")
    print("Skipping HWE filter")
    # paste0(round((table(sample_gt_renamed$HWE < 0.05)[[2]]/quilt_variants*100),2),"%")
  } else if(cross_type == "het3"){
    print("HET3 strains; sites not expected to adhere to HWE")
    print("Skipping HWE filter")
    # paste0(round((table(sample_gt_renamed$HWE < 0.05)[[2]]/quilt_variants*100),2),"%")
  } else {
    print("Pct of sites that deviate from HWE:")
    # paste0(round((table(sample_gt_renamed$HWE < 0.05)[[2]]/quilt_variants*100),2),"%")

    # filtering by HWE
    sample_gt_renamed <- sample_gt_renamed[which(sample_gt_renamed$HWE > 0.05),]
  }
} else {
  print("Chromosome X; sites not necessarily expected to adhere to HWE")
  print("Skipping HWE filter")
  # paste0(round((table(sample_gt_renamed$HWE < 0.05)[[2]]/quilt_variants*100),2),"%")
}

# filtering by info score
lower_info_score <- 0.95
above_threshold_sites <- length(which(sample_gt_renamed$INFO_SCORE > lower_info_score))
paste0(signif(above_threshold_sites/length(sample_gt_renamed$INFO_SCORE), 4)*100,"% of sites above 0.95 INFO score threshold (",above_threshold_sites,")")

# This bin is for runs where there were sites above the threshold, but fewer
message(paste0(above_threshold_sites," info scores above 0.95"))
sample_gt_renamed = sample_gt_renamed[which(sample_gt_renamed$INFO_SCORE > lower_info_score),]
filtered_quilt_variants <- nrow(sample_gt_renamed)
new_threshold <- lower_info_score

# apply the changes to founder_gt
founder_gt_renamed = founder_gt_renamed[rownames(founder_gt_renamed) %in% rownames(sample_gt_renamed),]
filtered_quilt_variants <- nrow(sample_gt_renamed)
resolution_summary <- data.frame(quilt_variants,filtered_quilt_variants, new_threshold)
resolution_summary$chr <- chr

# load 1M grid
grid <- read.csv(gridfile)
chr_grid <- grid[which(grid$chr == chr),]

# find nearest sample gt to grid locations
nearest_base <- function(query, subject) {
  # Ensure query and subject are data frames with columns 'start' and 'end'
  if (!all(c("start", "end") %in% colnames(query)) || !all(c("start", "end") %in% colnames(subject))) {
    stop("Both query and subject must have 'start' and 'end' columns")
  }
  
  nearest_indices <- sapply(query$start, function(q_start) {
    distances <- abs(subject$start - q_start)
    nearest_index <- which.min(distances)
    return(nearest_index)
  })
  
  return(nearest_indices)
}
sample_subject <- data.frame(start = sample_gt_renamed$POS, end = sample_gt_renamed$POS)
grid_query <- data.frame(start = chr_grid$pos*1e6, end = chr_grid$pos*1e6)
nearest_indices <- nearest_base(grid_query, sample_subject)
nearest_indices <- unique(nearest_indices)

# anchor to the grid
anchored_sample_gt <- sample_gt_renamed[nearest_indices,]
anchored_founder_gt <- founder_gt_renamed[nearest_indices,]
stopifnot(rownames(anchored_sample_gt) == rownames(anchored_founder_gt))
resolution_summary$grid_variants <- nrow(anchored_sample_gt)
write.csv(resolution_summary, paste0("chr",chr,"_resolution_summary.csv"), quote = F, row.names = F)

# Merge founders and samples together.
all_gt <- dplyr::full_join(anchored_founder_gt, anchored_sample_gt) %>%
  dplyr::select(-HWE, -INFO_SCORE) %>%
  dplyr::arrange(CHROM, POS)
rownames(all_gt) <- rownames(anchored_founder_gt)

message("Get allele codes for each SNP")
# Get the allele codes for each SNP.
ref = as.character(all_gt$REF)
alt = as.character(all_gt$ALT)
alleles = cbind(ref, alt)
rm(ref, alt)

#Encode the genotypes from qtl2.
all_gt_meta <- all_gt[,1:4]
if(cross_type == "bxd"){
  all_gt = encode_geno(all_gt[,-c(1:4)], alleles, output_codes = c("-","B","H","D"))
} else {
  all_gt = encode_geno(all_gt[,-c(1:4)], alleles)
}
all_gt <- cbind(all_gt_meta, all_gt)


# format founder genos for qtl2
founder_gt <- all_gt %>%
  dplyr::select(colnames(founder_gt_renamed)[-c(1:4)])

if(cross_type == "do"){
  
  colnames(founder_gt) = LETTERS[1:ncol(founder_gt)]
  founder_gt = data.frame(marker = rownames(founder_gt), founder_gt)
  
} else if(cross_type == "bxd"){
  
  colnames(founder_gt) = c("B","D")
  founder_gt = data.frame(marker = rownames(founder_gt), founder_gt)
  
} else if(cross_type == "het3"){
  
  colnames(founder_gt) = c("Y","B","C","D")
  founder_gt = data.frame(marker = rownames(founder_gt), founder_gt)
  
} else {
  
  colnames(founder_gt) = LETTERS[1:ncol(founder_gt)]
  founder_gt = data.frame(marker = rownames(founder_gt), founder_gt)
  
}
write.csv(founder_gt, file = paste0("chr",chr,"_founder_geno.csv"),
          quote = FALSE, row.names = FALSE)

# Write out sample genotypes.
sample_gt <- all_gt %>%
  dplyr::select(colnames(anchored_sample_gt)[-c(1:6)])
sample_gt = data.frame(marker = rownames(sample_gt),
                       sample_gt, check.names = F)
write.csv(sample_gt, file = paste0("chr",chr,"_sample_geno.csv"),
          quote = FALSE, row.names = FALSE)

# Write out physical map.
pmap = data.frame(marker = rownames(anchored_founder_gt),
                  chr    = anchored_founder_gt$CHROM,
                  pos    = anchored_founder_gt$POS * 1e-6)
write.csv(pmap, file = paste0("chr",chr,"_pmap.csv"),
          quote = FALSE, row.names = FALSE)

# Write out genetic map.
markers = read.delim(marker_file, sep = ' ')
markers = subset(markers, position %in% anchored_founder_gt$POS)
stopifnot(markers$position == anchored_founder_gt$POS)
gmap = data.frame(marker = rownames(anchored_founder_gt),
                  chr    = anchored_founder_gt$CHROM,
                  pos    = markers$Genetic_Map.cM.)
write.csv(gmap, file = paste0("chr",chr,"_gmap.csv"),
          quote = FALSE, row.names = FALSE)


# Write out covariates

# Has a cross type been specified?
stopifnot("cross_type" %in% ls(pattern = "cross_type"))

# Does the covar file match the format of the cross type?
if(cross_type == "genail4" | cross_type == "genail8" | cross_type == "cc" | cross_type == "het3" | cross_type == "F1_mut"){
  stopifnot("id" %in% colnames(meta))
  stopifnot("gen" %in% colnames(meta))
  stopifnot(colnames(founder_gt)[-1] %in% LETTERS)
  if("sex" %in% colnames(meta)){
    meta$sex[meta$sex == "female" | meta$sex == "f"] <- "F"
    meta$sex[meta$sex == "male" | meta$sex == "m"] <- "M"
  } else {
    message("No sex column included")
  }
  write.csv(meta, file = 'covar.csv',
            quote = FALSE, row.names = FALSE)

} else if(cross_type == "do"){
  stopifnot("id" %in% colnames(meta))
  stopifnot("gen" %in% colnames(meta))
  if("sex" %in% colnames(meta)){
    meta$sex[meta$sex == "female" | meta$sex == "f"] <- "F"
    meta$sex[meta$sex == "male" | meta$sex == "m"] <- "M"
  } else {
    message("No sex column included")
  }
  write.csv(meta, file = 'covar.csv',
            quote = FALSE, row.names = FALSE)

} else if(cross_type == "bxd"){
  
  # keep the default id column
  stopifnot("id" %in% colnames(meta))
  if("sex" %in% colnames(meta)){
    meta$sex[meta$sex == "female" | meta$sex == "f"] <- "F"
    meta$sex[meta$sex == "male" | meta$sex == "m"] <- "M"
  }
  write.csv(meta, file = 'covar.csv',
            quote = FALSE, row.names = FALSE)

} else {
  "Cross type parsing not included for the specified cross type"
}


# Write out phenotypes.
pheno = data.frame(id  = meta$id,
                   val = rep(1, nrow(meta)))
write.csv(pheno, file = 'pheno.csv',
           quote = FALSE, row.names = FALSE)

