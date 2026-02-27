################################################################################
# Given a counts file, use Ddx3y and Xist to determine the sex of each sample
# in the counts file. Generate a sample metadata file containing the sample ID
# and the sex of each mouse.
#
# Arguments:
# arg 1: string containing the full path to the counts file for one sample.
#        This should be a tab-delimited file output from the CS-NF RNASeq
#        pipeline. The first column should be named "gene_id" and 
#        should contain the Ensembl ID, gene symbol, or a concatenated
#        <Ensembl ID>_<gene symbol> name. There must also be a column
#        called "TPM".
# arg 2: string containing the full path to the directory to write files.
# arg 3: string containing the file name prefix for output file.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2024-09-23
################################################################################

##### VARIABLES #####

args = commandArgs(trailingOnly = TRUE)

if(length(args) != 3) {

  stop('ERROR: Three arguments required in determine.sex.R.\nusage: Rscript determine_sex.R <counts_file> <output_path> <file_prefix>')

} # if(length(args) != 3)

# Path to the counts file.
counts_file = args[1]

# Path to the output directory.
output_path = args[2]

# Prefix for output file.
file_prefix = args[3]

# Mouse and human Ensembl IDs & symbols.
## If an additional annotation set is added, add the relevant IDs to the grep strings.
xist_search_str  = 'ENSMUSG00000086503|ENSG00000229807|xist|MGI:98974|Xist'
ddx3y_search_str = 'ENSMUSG00000069045|ENSG00000067048|ddx3y|MGI:1349406|Ddx3y'

# Note on these genes: 
## Ddx3y is a Y-linked gene that should only be expressed in males. 
## Xist is an X-linked gene that is typically expressed at higher levels in females due to X-chromosome inactivation. 
## By comparing the expression levels of these two genes, we can infer the sex of the sample. 
## If Ddx3y is expressed and Xist is not, the sample is likely male. 
## If Xist is expressed and Ddx3y is not, the sample is likely female.
## If both genes are expressed at similar levels, the sample may be ambiguous or of low quality, and further investigation may be needed.

sex_threshold <- 0.1
# This threshold is used to determine if the expression levels of Xist and Ddx3y are close enough to be considered "undetermined".

##### FUNCTIONS #####

# Given an Ensembl gene ID and a gene symbol, search for the row which
# contains one or the other in column 1 of the counts. 
# Arguments:
# search_str: string containing the Ensembl ID and/or gene symbol of of gene 
#             to search for. e.g. "ENSMUSG00000086503|ENSG00000229807|xist"
# counts:  data.frame containing the gene counts. First column MUST
#          contain either the Ensembl ID or gene symbols 
# Returns: integer that is the row containing the requested gene.
find_gene_row = function(search_str, counts) {

  message('[INFO] Searching for pattern: ', search_str)

  # Search for the gene.
  gene_row = grep(search_str, counts$gene_id, ignore.case = TRUE)

  # If we didn't find it, quit with an error.
  # TBD: Return -1 and let the caller handle errors?
  if(length(gene_row) == 0) {

    message('[ERROR] No match found for pattern: ', search_str)
    error_msg = paste('determine_sex.R: Could not find', search_str, 
                      'in counts file.')
    message(error_msg)
    
    # Write output file with NA values
    out_file = file.path(output_path, paste0(file_prefix, '_SexDet.csv'))
    cat(error_msg, file = out_file, sep = "\n")
    message('[INFO] Output written to: ', out_file, ' (no gene match found)')
    
    quit(status = 0, save = "no")
  
  } # if(length(gene_row) == 0)

  # Warn if more than one row matched — likely an annotation collision.
  if(length(gene_row) > 1) {
    message('[WARN] Multiple rows matched pattern "', search_str, '": rows ',
            paste(gene_row, collapse = ', '), ' — using first match only.')
    message('[WARN] Matched gene_id values: ',
            paste(counts$gene_id[gene_row], collapse = ' | '))
    gene_row = gene_row[1]
  } # if(length(gene_row) > 1)

  message('[INFO] Match found at row ', gene_row,
          ' | gene_id: ', counts$gene_id[gene_row])

  return(gene_row)

} # find_gene_row()


##### MAIN #####

# Verify that the file exists.
if(!file.exists(counts_file)) {

  stop(paste('ERROR: File', counts_file, 'not found.'))

} # if(!file.exists(counts_file))

message('[INFO] Reading counts file: ', counts_file)

# Read in the counts.
counts = read.delim(counts_file)

message('[INFO] Loaded ', nrow(counts), ' rows x ', ncol(counts), ' columns from counts file.')

# Verify that we have a "gene_id" column.
if(!'gene_id' %in% colnames(counts)) {

  stop(paste('determine_sex.R: column name gene_id not found in', counts_file))

} # if(!'gene_id' %in% colnames(counts))

# Verify that we have a "TPM" column.
if(!'TPM' %in% colnames(counts)) {

  stop(paste('determine_sex.R: column name TPM not found in', counts_file))

} # if(!'TPM' %in% colnames(counts))

# Find the row conatining Xist counts.
xist_row = find_gene_row(search_str = xist_search_str,
                         counts     = counts)

# Find the row conatining Ddx3y counts.
ddx3y_row = find_gene_row(search_str = ddx3y_search_str,
                          counts     = counts)

# Get the counts for the two genes.
sex_counts        = counts$TPM[c(xist_row, ddx3y_row)]
names(sex_counts) = c('Xist', 'Ddx3y')

message('[INFO] TPM values — Xist: ', sex_counts['Xist'], ' | Ddx3y: ', sex_counts['Ddx3y'])

sex = data.frame(id    = file_prefix,
                 xist  = sex_counts['Xist'],
                 ddx3y = sex_counts['Ddx3y'],
                 sex   = ifelse(sex_counts['Xist'] == 0 & sex_counts['Ddx3y'] < sex_threshold,
                                NA, 
                                ifelse(abs(sex_counts['Xist'] - sex_counts['Ddx3y']) / 
                                       (sex_counts['Xist'] + sex_counts['Ddx3y'] + 1) < 0.1,
                                       'Undetermined',
                                       ifelse(sex_counts['Xist'] > sex_counts['Ddx3y'], 
                                              'female', 'male'))
                                )
) 

message('[INFO] Sex determination result — sample: ', file_prefix,
        ' | xist_TPM: ', sex$xist,
        ' | ddx3y_TPM: ', sex$ddx3y,
        ' | sex: ', sex$sex)

# Write out a file containing the estimated sex for each sample.
out_file = file.path(output_path, paste0(file_prefix, '_SexDet.csv'))
write.csv(sex, file = out_file, quote = FALSE, row.names = FALSE)
message('[INFO] Output written to: ', out_file)
