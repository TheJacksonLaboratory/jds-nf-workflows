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
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)



## set command line arguments ----
args <- commandArgs(trailingOnly = TRUE)

# Wrapped fasta file : Mus_musculus.GRCm38.cds.all.LINEWRAP.fa
wrap_fa = args[1]
# gen_org = args[2]

suffix   = strsplit(wrap_fa, '.', fixed = TRUE)[[1]][1]
suffix_2 = strsplit(wrap_fa, '.', fixed = TRUE)[[1]][2]

# Scan through the fasta file to get transcript names and lengths
transcripts <- scanFasta(wrap_fa)
nsequences <- nrow(transcripts) - sum(transcripts$Duplicate)


# Assign a random TPM value to each non-duplicated transcript sequence
TPMs <- rep(0, nrow(transcripts))
TPMs[!transcripts$Duplicate] <- rexp(nsequences) # those that are duplicated are given '0' TPM values. 


# PE
lib_size = c(1000000, 30000000) 

for(i in seq_along(lib_size)) {

  print(lib_size[i])

  if(lib_size[i] == 1000000){
    reads = '1m'
  } else
    reads = '30m'

  print(reads)

  output   = paste0(suffix, '_simRNA_', reads,'Reads') 
  output_2 = paste0(suffix, '.', suffix_2, '.', reads,'_isoform_true_expression.txt') 

  print(output)
  print(output_2)
  print(wrap_fa)
  print(lib_size[i])  

  true.counts <- simReads(

    # the transcript database and the wanted abundances
    wrap_fa,TPMs, output, 
    # options on the output
    library.size = lib_size[i],
    read.length = 100,
    truth.in.read.names = FALSE,
  
    # simulating sequencing errors
    simulate.sequencing.error = TRUE,
    quality.reference = NULL,
  
    # parameters for generating paired-end reads.
    paired.end = TRUE,
    fragment.length.min = 125,
    fragment.length.max = 500,
    fragment.length.mean = 250,
    fragment.length.sd = 40,
  
    # manipulating transcript names
    simplify.transcript.names = FALSE)

write.table(x = true.counts, file = output_2, sep = '\t', quote = F, row.names = F)

}


# SE for reads 1m
reads = '1m'

output   = paste0(suffix, '_simRNA_', reads,'Reads_SE') 
output_2 = paste0(suffix, '.', suffix_2, '.', reads,'_isoform_SE_true_expression.txt') 

print(output)
print(output_2)
print(wrap_fa)

true.counts <- simReads(

    # the transcript database and the wanted abundances
    wrap_fa,TPMs, output,
    # options on the output
    library.size = lib_size[1],
    read.length = 100,
    truth.in.read.names = FALSE,

    # simulating sequencing errors
    simulate.sequencing.error = TRUE,
    quality.reference = NULL,

    # parameters for generating paired-end reads.
    paired.end = FALSE,
    fragment.length.min = 125,
    fragment.length.max = 500,
    fragment.length.mean = 250,
    fragment.length.sd = 40,

    # manipulating transcript names
    simplify.transcript.names = FALSE)

write.table(x = true.counts, file = output_2, sep = '\t', quote = F, row.names = F)





