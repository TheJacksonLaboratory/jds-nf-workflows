## BJS Note: this script was located in the root path of the Docker container
## gcr.io/nygc-public/sv_cnv@sha256:1c14a50d131323a2a4bab323cf224879776af8de37f93df79292fd2e63269274
## It is reproduced below as it exists there without modification

## Filter a bedpe for somatic variants (i.e., not in specified germline databases), and 
## high-confidence variants (2+ callers or 1 caller with a nearby copy number changepoint)
libs = c('optparse', 'GenomicRanges', 'dplyr', 'tidyverse', 'stringr')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T, quietly=T)))
options(width=200, scipen=999)

## Check if databases db are in info string x
inDatabase = function(x, db) {
  
  ## Split info field, look for database entry
  x = unlist(strsplit(x, ';', fixed=T))
  x = grep('known=', x, fixed=T, value=T)
  x = gsub('known=', '', x, fixed=T)
  x = unlist(strsplit(x, ',', fixed=T))
  
  return(any(x %in% db))
  
}

## Check if databases db are in info string x
inDatabaseMatchType = function(x, db) {

  x$ID <- seq.int(nrow(x))

  x <- left_join(x, x %>%
  separate_rows(info, sep = ";") %>%
  filter(str_detect(info, "known=")) %>%
  mutate(info = gsub('known=', '', info),
         to_filter = ifelse(str_detect(info, paste(!!db, collapse = "|")) & str_detect(info, type), 1, 0 )) %>% # collapse converts a pasted array to | delim string needed as 'or' in str_detect
  dplyr::select(ID, to_filter),
  by = 'ID') %>%
  dplyr::filter(to_filter != 1) %>%
  dplyr::select(-ID, -to_filter)

  return(x)
  
}

makeGRangesFromChangepoint = function(x) {
  
  x = unlist(strsplit(x, ':|-'))
  GRanges(seqnames=x[1], ranges=IRanges(as.numeric(x[2:3]), as.numeric(x[2:3])))
                                              
}

## Is variant x a high-confidence variant? 
## Meant to be used with apply(,2,)
isHighConfidence = function(x) {
  
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
  make_option(c("-f", "--filter_databases"),          type='character',  help="Comma-separated list of databases to filter, looking in the info field"),
  make_option(c("-s", "--out_file_somatic"),          type='character',  help="Output somatic BEDPE"),
  make_option(c("-o", "--out_file_highconf"),         type='character',  help="Output high-confidence BEDPE"),
  make_option(c("-g", "--genome"),                    type='character',  help="Genome version", default = "GRCh38"))
opt = parse_args(OptionParser(option_list=option_list))

## Unpack arguments
opt$filter_databases = unlist(strsplit(opt$filter_databases, ',', fixed=T))

## Read bedpe
x = read.csv(opt$bedpe, h=T, stringsAsFactors=F, sep='\t', check.names=F)

## Filter for known germline variants 
if (opt$genome == 'GRCh38') {
  x = x[!sapply(x$info, inDatabase, opt$filter_databases), ]
} else if (opt$genome == 'GRCm39') {
  x = inDatabaseMatchType(x, opt$filter_databases)
} # for GRCm39 we know the type of SV event, and only wish to filter cases where the SV event matches the database hit. 

## Write out somatic variants
write.table(x, opt$out_file_somatic, row.names=F, col.names=T, sep='\t', quote=F)

## Filter for high confidence
x = x[apply(x, 1, isHighConfidence), ]

## Write result
write.table(x, opt$out_file_highconf, row.names=F, col.names=T, sep='\t', quote=F)
