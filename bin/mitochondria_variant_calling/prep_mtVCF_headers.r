# Required libraries
suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(dplyr)
    library(optparse)
})

option_list <- list(
    make_option(c("-v", "--vcf"), type = "character", help = "VCF files"),
    make_option(c("-c", "--caller"), type = "character", help = "Caller name for VCF INFO field"),
    make_option(c("-s", "--sampleID"), type = "character", help = "Sample ID for reheader"),
    make_option(c("-o", "--output"), type = "character", default = "reheaded-file.vcf", help = "Output VCF file name [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$vcf) || is.null(opt$caller) || is.null(opt$sampleID)) {
    stop("--vcf, --caller, and --sampleID options must be provided.")
}

vcf_files_vec <- strsplit(opt$vcf, ",")[[1]]
caller_names_vec <- strsplit(opt$caller, ",")[[1]]

if (length(vcf_files_vec) != length(caller_names_vec)) {
    stop("Number of VCF files and caller names must match.")
}
if (length(vcf_files_vec) != 1) {
    stop("Exactly one VCF file must be provided.")
}
if (!file.exists(vcf_files_vec)) {
    stop(paste("VCF file does not exist:", vcf_files_vec))
}

# Read VCF
vcf <- readVcf(opt$vcf)

# Add a new INFO field to the header if not present
hdr <- header(vcf)

# Add INFO fields for CALLER and SUPPORT if not present
if (!"CALLER" %in% rownames(info(hdr))) {
    info(hdr) <- rbind(
        info(hdr),
        DataFrame(
            Number = "1",
            Type = "String",
            Description = "Variant caller source",
            row.names = "CALLER"
        )
    )
    header(vcf) <- hdr
}
if (!"SUPPORT" %in% rownames(info(hdr))) {
    info(hdr) <- rbind(
        info(hdr),
        DataFrame(
            Number = "1",
            Type = "Integer",
            Description = "Variant support evidence",
            row.names = "SUPPORT"
        )
    )
    header(vcf) <- hdr
}

# Ensure FORMAT$AF is Number=A and Type=Float
if ("AF" %in% rownames(geno(hdr))) {
    if (geno(hdr)["AF", "Number"] != "A") {
        geno(hdr)["AF", "Number"] <- "A"
    }
    if (geno(hdr)["AF", "Type"] != "Float") {
        geno(hdr)["AF", "Type"] <- "Float"
    }
    header(vcf) <- hdr
}

# Add the caller name to the INFO field
info(vcf)$CALLER <- opt$caller
info(vcf)$SUPPORT <- 1
# Change genotype/sample name to opt$sampleID
colnames(vcf) <- opt$sampleID

suppressWarnings(writeVcf(vcf, filename = opt$output, bgzip = TRUE, index = TRUE))
