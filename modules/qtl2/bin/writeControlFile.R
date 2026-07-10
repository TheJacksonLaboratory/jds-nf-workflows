#!/usr/bin/env Rscript
library(qtl2)
library(qtl2convert)
library(dplyr)
################################################################################
# Make .json file and cross object from genotype files and GigaMUGA reference
# data.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20250521
################################################################################ 

# testing
# test_dir = "/flashscratch/widmas/HR_QC_outputDir/work/1c/5ccd8e7e35c32b08e88e6f5a4a1fb4/"
# setwd(test_dir)

# covar file
covar_file <- list.files(pattern = "covar.csv")

cat(" --Current working directory:")
print(getwd())

# read in covariate file and make sex codes uniform
covar <- read.csv(covar_file, stringsAsFactors = F)
covar$sex[covar$sex == "female" | covar$sex == "f" | covar$sex == "FALSE" | covar$sex == FALSE | covar$sex == "XX" | covar$sex == "XO"] <- "F"
covar$sex[covar$sex == "male" | covar$sex == "m" | covar$sex == "XY" | covar$sex == "XXY"] <- "M"

# write out new covar file
covar <- covar %>%
  dplyr::select(id, sex, gen, everything()) %>%
  dplyr::distinct(id, sex, gen)
write.csv(x = covar, file = "sex_corrected_covar.csv", row.names = F, quote = F)

# move elements to one directory for reference and for proper cross paths
ostem = file.path(getwd(),"out")
dir.create(ostem, showWarnings = T, recursive = T)
system(paste0("mv GM*.csv ",ostem))
system(paste0("mv geno*.csv ",ostem))
system(paste0("mv sex_corrected_covar.csv ",ostem))

# Write control file
chr <- c(1:19, "X")
cat(" -Writing control file\n")
# for now stick with "do" cross type, eventually take a crosstype param
if(length(unique(covar$sex)) > 1){
  qtl2::write_control_file(output_file = file.path(ostem,"QC_HAP.json"),
                           crosstype="do",
                           founder_geno_file=paste0("GM_foundergeno", chr, ".csv"),
                           founder_geno_transposed=TRUE,
                           gmap_file=paste0("GM_gmap", chr, ".csv"),
                           pmap_file=paste0("GM_pmap", chr, ".csv"),
                           geno_file=paste0("geno", chr, ".csv"),
                           geno_transposed=TRUE,
                           geno_codes=list(A=1, H=2, B=3),
                           xchr="X",
                           covar_file="sex_corrected_covar.csv",
                           crossinfo_covar="gen",
                           sex_covar="sex",
                           sex_codes=c("F"="female", "M"="male"),
                           overwrite = T)
} else if(unique(covar$sex) == "F"){
  cat(" Note: only females present in cross.\n")
  qtl2::write_control_file(output_file = file.path(ostem,"QC_HAP.json"),
                           crosstype="do",
                           # description="QC_HAP",
                           founder_geno_file=paste0("GM_foundergeno", chr, ".csv"),
                           founder_geno_transposed=TRUE,
                           gmap_file=paste0("GM_gmap", chr, ".csv"),
                           pmap_file=paste0("GM_pmap", chr, ".csv"),
                           geno_file=paste0("geno", chr, ".csv"),
                           geno_transposed=TRUE,
                           geno_codes=list(A=1, H=2, B=3),
                           xchr="X",
                           covar_file="sex_corrected_covar.csv",
                           crossinfo_covar="gen",
                           sex_covar="sex",
                           sex_codes=c("F"="female"),
                           overwrite = T)
} else {
  cat(" Note: only males present in cross.\n")
  qtl2::write_control_file(output_file = file.path(ostem,"QC_HAP.json"),
                           crosstype="do",
                           # description="QC_HAP",
                           founder_geno_file=paste0("GM_foundergeno", chr, ".csv"),
                           founder_geno_transposed=TRUE,
                           gmap_file=paste0("GM_gmap", chr, ".csv"),
                           pmap_file=paste0("GM_pmap", chr, ".csv"),
                           geno_file=paste0("geno", chr, ".csv"),
                           geno_transposed=TRUE,
                           geno_codes=list(A=1, H=2, B=3),
                           xchr="X",
                           covar_file="sex_corrected_covar.csv",
                           crossinfo_covar="gen",
                           sex_covar="sex",
                           sex_codes=c("M"="male"),
                           overwrite = T)
}

 
cat(" -Reading control file\n")
# read cross
cross <- qtl2::read_cross2(file = file.path(ostem,"QC_HAP.json"), quiet = F)

# sample hash to link to intensity input
int_dedup_file <- "dedup_samples.csv"
hash <- strsplit(covar_file,"_")[[1]][[1]]
dedup_samples <- read.csv(int_dedup_file)
dedup_samples <- dedup_samples$sample[which(dedup_samples$hash == hash)]

# subset cross to individuals from the correct neogen file if sample was genotyped multiple times
cross <- subset(cross, ind = dedup_samples)

saveRDS(cross, file = "preQC_cross.rds")

