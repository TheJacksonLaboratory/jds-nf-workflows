################################################################################
# 
#
# simReads: Generate Simulated Reads from a Set of Transcripts 
#
# The simReads function generates simulated reads from a set of transcripts. 
# Each transcript has a predefined abundance in the output.
#
#
# Container: simreads.sif
#
#
# Sai Lek
# sai.lek@jax.org
# 2026-03-11
################################################################################

##### LIBRARIES #####

library(Rsubread)

## set command line arguments ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 10) {
  stop(
    "Usage: Rscript generate_simulated_RNA_data.R <wrapped_fasta.fa> <library_size[,library_size...]> <PE|SE|BOTH> <read_length> <fragment_length_min> <fragment_length_max> <fragment_length_mean> <fragment_length_sd> <simulate_sequencing_error> <quality_reference>",
    call. = FALSE
  )
}

# Wrapped fasta file : Mus_musculus.GRCm38.cds.all.LINEWRAP.fa
wrap_fa <- args[1]

# Accept one value (e.g. 1000000) or comma-separated values (e.g. 1000000,30000000)
library_size <- as.numeric(trimws(strsplit(args[2], ",", fixed = TRUE)[[1]]))
if (any(is.na(library_size)) || any(library_size <= 0)) {
  stop("library_size must be a positive integer or comma-separated positive integers.", call. = FALSE)
}

library_strategy <- toupper(trimws(args[3]))
if (!(library_strategy %in% c("PE", "SE", "BOTH"))) {
  stop("library_strategy must be one of: PE, SE, BOTH", call. = FALSE)
}

read_length <- as.numeric(args[4])
fragment_length_min <- as.numeric(args[5])
fragment_length_max <- as.numeric(args[6])
fragment_length_mean <- as.numeric(args[7])
fragment_length_sd <- as.numeric(args[8])
simulate_sequencing_error <- as.logical(args[9])
quality_reference <- args[10]


suffix <- strsplit(wrap_fa, ".", fixed = TRUE)[[1]][1]
suffix_2 <- strsplit(wrap_fa, ".", fixed = TRUE)[[1]][2]

# Scan through the fasta file to get transcript names and lengths
transcripts <- scanFasta(wrap_fa)
nsequences <- nrow(transcripts) - sum(transcripts$Duplicate)

# Assign a random TPM value to each non-duplicated transcript sequence
TPMs <- rep(0, nrow(transcripts))
TPMs[!transcripts$Duplicate] <- rexp(nsequences) # those that are duplicated are given '0' TPM values. 

read_label <- function(size_value) {
  if (size_value %% 1000000 == 0) {
    paste0(size_value / 1000000, "m")
  } else {
    as.character(size_value)
  }
}

run_simulation <- function(size_value, paired_end_value) {
  reads <- read_label(size_value)

  if (paired_end_value) {
    output <- paste0(suffix, "_simRNA_", reads, "Reads")
    output_2 <- paste0(suffix, ".", suffix_2, ".", reads, "_isoform_true_expression.txt")
  } else {
    output <- paste0(suffix, "_simRNA_", reads, "Reads_SE")
    output_2 <- paste0(suffix, ".", suffix_2, ".", reads, "_isoform_SE_true_expression.txt")
  }

  true.counts <- simReads(

    # the transcript database and the wanted abundances
    wrap_fa, TPMs, output,
    # options on the output
    library.size = size_value,
    read.length = read_length,
    truth.in.read.names = FALSE,

    # simulating sequencing errors
    simulate.sequencing.error = simulate_sequencing_error,
    quality.reference = if (simulate_sequencing_error && !quality_reference %in% c("none", "", NA)) quality_reference else NULL,

    # parameters for generating paired-end reads.
    paired.end = paired_end_value,
    fragment.length.min = fragment_length_min,
    fragment.length.max = fragment_length_max,
    fragment.length.mean = fragment_length_mean,
    fragment.length.sd = fragment_length_sd,

    # manipulating transcript names
    simplify.transcript.names = FALSE
  )

  write.table(x = true.counts, file = output_2, sep = "\t", quote = FALSE, row.names = FALSE)
}

for (size_value in library_size) {
  if (library_strategy %in% c("PE", "BOTH")) {
    run_simulation(size_value = size_value, paired_end_value = TRUE)
  }
  if (library_strategy %in% c("SE", "BOTH")) {
    run_simulation(size_value = size_value, paired_end_value = FALSE)
  }
}
