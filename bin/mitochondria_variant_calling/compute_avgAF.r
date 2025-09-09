# Required libraries
suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(dplyr)
    library(optparse)
})

option_list <- list(
    make_option(c("-v", "--vcf"), type = "character", help = "VCF files"),
    make_option(c("-s", "--sampleID"), type = "character", help = "Sample ID for reheader"),
    make_option(c("-o", "--output"), type = "character", default = "reheaded-file.vcf", help = "Output VCF file name [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$vcf) || is.null(opt$sampleID)) {
    stop("--vcf and --sampleID options must be provided.")
}

vcf_files_vec <- strsplit(opt$vcf, ",")[[1]]

if (length(vcf_files_vec) != 1) {
    stop("Exactly one VCF file must be provided.")
}   

if (!file.exists(vcf_files_vec)) {
    stop(paste("VCF file does not exist:", vcf_files_vec))
}

# Read VCF
vcf <- readVcf(opt$vcf, row.names = FALSE)

# Add a new INFO field to the header if not present
hdr <- header(vcf)

# Add INFO fields for CALLER and SUPPORT if not present
if (!"AVG_AF" %in% rownames(info(hdr))) {
    info(hdr) <- rbind(
        info(hdr),
        DataFrame(
            Number = "A",
            Type = "Float",
            Description = "Average allele frequency across callers for a particular variant",
            row.names = "AVG_AF"
        )
    )
    header(vcf) <- hdr
}

# Compute average allele frequency (AVG_AF) for each variant
af_matrix <- geno(vcf)$AF

# If AF is missing, set AVG_AF to NA
if (is.null(af_matrix)) {
    avg_af <- rep(NA_real_, nrow(vcf))
} else {
    # AF can be a matrix (variants x samples or alleles)
    # We want to average across all available AF values for each variant
    avg_af <- apply(af_matrix, 1, function(x) {
        # x may be a vector (per sample or allele), average non-NA values
        if (all(is.na(x))) {
            NA_real_
        } else {
            round(mean(as.numeric(x), na.rm = TRUE), 4)
        }
    })
}

# Set the INFO field AVG_AF
info(vcf)$AVG_AF <- avg_af

writeVcf(vcf, filename = opt$output)
