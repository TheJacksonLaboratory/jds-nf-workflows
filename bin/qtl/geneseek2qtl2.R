#!/usr/bin/env Rscript
# load required packages
library(data.table)
library(qtl2convert)
library(vroom)
library(parallel)
library(dplyr)
library(fst)
library(purrr)
args <- commandArgs(trailingOnly = TRUE)
cat(paste0(" -- R/qtl2 version: ",qtl2::qtl2version(),"\n"))
cat(paste0(" -- R/qtl2convert version: ",packageVersion("qtl2convert"),"\n"))

# allele codes
codefile <- args[1]
# codefile <- "/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/bin/CC_DO_data/GM_allelecodes.csv"


# metadata file indicating which mice should be retained from array files
metadata_path <- args[2]
# metadata_path <- "/projects/compsci/vmp/USERS/widmas/haplotype_reconstruction_qtl-nf/sample_sheets/20250509_attie_covar.csv"

# FinalReport files
ifile <- args[3]
# ifile <- "/projects/compsci/vmp/USERS/widmas/attie_500/data/genotypes/Univ_of_Wisconsin_Schueler_MURGIGV01_20221021/Univ_of_Wisconsin_Schueler_MURGIGV01_20221021_FinalReport.zip"

cat("Reading covar file")
metadata <- read.csv(metadata_path, tryLogical = F)
metadata$id <- as.character(metadata$id)
if("include" %in% colnames(metadata)){
  cat(" -'include' field detected in metadata file.\n")
  cat(" -Filtering metadata to included samples\n")
  metadata <- metadata %>%
    dplyr::filter(include == TRUE)
}

# read genotype codes
codes <- data.table::fread(codefile, skip = 3, data.table = F)

# # make excluded samples object
# excluded_samples <- NULL

## MAIN ##
cat(" -Reading data\n")
# read in the genotypes
g <- suppressMessages(vroom::vroom(file = ifile, 
                                   skip = 9, 
                                   num_threads = parallel::detectCores()))
  
# make sure neogen samples are characters
if(!is.character(g$`Sample ID`)){
  g$`Sample ID` <- as.character(g$`Sample ID`)
}
  
# check that any of the samples in neogen file are in metadata
if(!any(unique(g$`Sample ID`) %in% metadata$id)){
  stop("No sample identifiers from supplied metadata match supplied FinalReport file(s). \nPlease sync file names. \n")
}
  
# filter neogen samples to those present in metadata
cat(" -Filtering genotypes to present samples\n")
g_2 <- g %>%
  dplyr::filter(`SNP Name` %in% codes$marker,
                `Sample ID` %in% metadata$id)
  
# NOTE: may need to revise the IDs in the 2nd column
samples <- unique(g_2$`Sample ID`)

# make matrix to contain the genotypes
geno <- matrix(nrow=nrow(codes), ncol=length(samples))
dimnames(geno) <- list(codes$marker, samples)

# fill in matrix
for(i in seq(along=samples)) {
  if(i==round(i,-1)) cat(" --Sample", i, "of", length(samples), "\n")
  wh <- (g_2$`Sample ID`==samples[i])
  geno[g_2[wh,]$`SNP Name`,i] <- paste0(g_2[wh,]$`Allele1 - Forward`, 
                                        g_2[wh,]$`Allele2 - Forward`)
}

cat(" -Encode genotypes\n")
geno <- qtl2convert::encode_geno(geno, as.matrix(codes[,c("A","B")]))
  
# make matrix of X intensities for sex checks
cat(" -Grab X and Y intensities\n")
gX <- g_2 %>%
  dplyr::filter(`SNP Name` %in% codes[which(codes$chr == "X"),]$marker,
                `Sample ID` %in% colnames(geno))
cX <- matrix(nrow=sum(codes$chr=="X"),
             ncol=length(samples))
dimnames(cX) <- list(codes[which(codes$chr == "X"),]$marker, samples)
  
# make matrix of Y intensities for sex checks
gY <- g_2 %>%
  dplyr::filter(`SNP Name` %in% codes[which(codes$chr == "Y"),]$marker,
                `Sample ID` %in% colnames(geno))
cY <- matrix(nrow=sum(codes$chr=="Y"),
             ncol=length(samples))
dimnames(cY) <- list(codes[which(codes$chr == "Y"),]$marker, samples)

for(i in seq(along=samples)) {
  if(i==round(i,-1)) cat(" --Sample", i, "of", length(samples), "\n")
  wh <- (gX[,"Sample ID"]==samples[i])
  cX[gX[wh,]$`SNP Name`,i] <- (gX$X[wh] + gX$Y[wh])/2
  
  wh <- (gY[,"Sample ID"]==samples[i])
  cY[gY[wh,]$`SNP Name`,i] <- (gY$X[wh] + gY$Y[wh])/2
}
  
# grab all marker intensities
cat(" -Grab all marker intensities\n")
g_3 <- g_2 %>%
  dplyr::select(`SNP Name`,`Sample ID`,X, Y) %>%
  dplyr::filter(`Sample ID` %in% colnames(geno))
colnames(g_3) <- c("snp","sample","X","Y")
int <- g_3 %>%
  tidyr::pivot_longer(-c(snp,sample), names_to = "channel") %>%
  tidyr::pivot_wider(names_from = sample, values_from = value) %>%
  data.frame(., check.names = F)


# check geno samples are intensity samples
stopifnot(all(colnames(int)[-c(1,2)] == colnames(geno)))
stopifnot(all(colnames(cX) == colnames(cY)))

# write X and Y intensities
cat(" -Writing X and Y intensities\n")
qtl2convert::write2csv(df = cbind(marker=rownames(cX), cX), 
                       filename = "chrXint.csv", 
                       comment = "X chr intensities",
                       overwrite=TRUE)
qtl2convert::write2csv(df = cbind(marker=rownames(cY), cY),
                       filename = "chrYint.csv",
                       comment ="Y chr intensities",
                       overwrite=TRUE)

# write data to chromosome-specific files
cat(" -Writing genotypes\n")
for(chrom in c(1:19,"X","Y","M")) {
  cat(" --Chromosome",chrom," \n")
  mar <- codes %>%
    dplyr::filter(chr == chrom) %>%
    dplyr::select(marker)
  g <- geno[rownames(geno) %in% mar$marker,]
  if(is.null(dim(g))){
    g <- as.matrix(g)
    colnames(g) <- colnames(geno)
  }
  qtl2convert::write2csv(df = cbind(marker=rownames(g), g),
                         filename = paste0("geno", chrom, ".csv"),
                         comment = paste0("genotypes for chr ", chrom),
                         overwrite=TRUE)
}

# write to fst file, maximally compressed
write_fst(int, "intensities.fst", compress=100)

# write the covariate file with a generic name
cat(" -Writing subsetted covar file \n")
filtered_meta <- metadata %>%
  dplyr::filter(id %in% colnames(geno)) %>%
  dplyr::mutate(provided_sex = sex)
write.csv(filtered_meta, "covar.csv", quote = F, row.names = F)
