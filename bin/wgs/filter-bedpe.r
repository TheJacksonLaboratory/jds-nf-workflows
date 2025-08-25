## BJS Note: this script was located in the root path of the Docker container
## gcr.io/nygc-public/sv_cnv@sha256:1c14a50d131323a2a4bab323cf224879776af8de37f93df79292fd2e63269274
## It is reproduced below as it exists there without modification

## Filter a bedpe for somatic variants (i.e., not in specified germline databases), and 
## high-confidence variants (2+ callers or 1 caller with a nearby copy number changepoint)
libs = c('optparse', 'GenomicRanges', 'dplyr', 'tidyverse', 'stringr')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T, quietly=T)))
options(width=200, scipen=999)

## Is variant x a high-confidence variant? 
## Meant to be used with apply(,2,)
isHighConfidence = function(x, cpmax) {
  
  ## Is there support from multiple callers?
  multi.caller = grepl(',', x['tools'])

  ## Is either breakpoint close enough to a changepoint?
  ## Changepoints are only annotated in the prior script if the 'nearest' neighbor is within a distance limit.
  ## Check if one or the other (or both) CNV change points exist and not NA or empty string
  if (is.null(x['cnv_changepoint_1']) || is.na(x['cnv_changepoint_1']) || x['cnv_changepoint_1'] == '') {
    near.ch1 = FALSE
  } else {
    near.ch1 = TRUE
  }

  if (is.null(x['cnv_changepoint_2']) || is.na(x['cnv_changepoint_2']) || x['cnv_changepoint_2'] == '') {
    near.ch2 = FALSE
  } else {
    near.ch2 = TRUE
  }
  
  return(multi.caller || near.ch1 || near.ch2)
  
}

## Collect arguments
option_list = list(
  make_option(c("-b", "--bedpe"),                     type='character',  help="Input BEDPE"),
  make_option(c("-o", "--outfile_highconf"),         type='character',  help="Output high-confidence BEDPE"))
opt = parse_args(OptionParser(option_list=option_list))

## Read bedpe
x = read.csv(opt$bedpe, h=T, stringsAsFactors=F, sep='\t', check.names=F)

## Filter for high confidence
x = x[apply(x, 1, isHighConfidence, opt$max_changepoint_distance), ]

## Write result
write.table(x, opt$outfile_highconf, row.names=F, col.names=T, sep='\t', quote=F)
