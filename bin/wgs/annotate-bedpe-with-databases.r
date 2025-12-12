## BJS Note: this script was located in the root path of the Docker container
## gcr.io/nygc-public/sv_cnv@sha256:1c14a50d131323a2a4bab323cf224879776af8de37f93df79292fd2e63269274
## It is reproduced below as it exists there without modification

## Annotate a merged bedpe with arbitrary databases
libs = c('optparse', 'StructuralVariantAnnotation', 'VariantAnnotation', 'rtracklayer', 'stringr', 'gUtils')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T, quietly=T)))
options(width=200, scipen=999)



## Handle non-standard bedpe columns better
readBEDPE = function(f) {
  
  ## Read file as Pairs object
  x = rtracklayer::import(f, format='bedpe')
  
  ## Update metadata column names
  x.mcol.names = colnames(read.csv(f, h=T, stringsAsFactors=F, sep='\t', check.names=F))
  colnames(mcols(x))[3:ncol(mcols(x))] = x.mcol.names[11:length(x.mcol.names)]
  mcols(x)$type = mcols(x)$name
  
  
  ## Convert to breakpoint ranges
  x  = StructuralVariantAnnotation::pairs2breakpointgr(x)
  
  return(x)
  
}



## Handle both BEDPE and BED files (headerless)
readDB = function(f) {
  
  is.bed = grepl('\\.bed\\.gz$|\\.bed$', f)
  is.bedpe = grepl('\\.bedpe\\.gz$|\\.bedpe$', f)
  is.vcf = grepl('\\.vcf\\.gz$|\\.vcf$', f)
  
  if (is.bed || is.bedpe) {
    x = rtracklayer::import(f)
    
    ## If this is a BEDPE, convert to a breakpointRanges object
    ## The package seems to misname SV type so update that
    if (is.bedpe) {
      x  = StructuralVariantAnnotation::pairs2breakpointgr(x)
      x$type = x$sourceId
    }
    
  }  
  
  if (is.vcf) {
    
    x = VariantAnnotation::readVcf(f)
        
    x = StructuralVariantAnnotation::breakpointRanges(x, nominalPosition=T)

    if ('svtype' %in% colnames(mcols(x))) {
      x$type = x$svtype
    } else {
      x$type = x$sourceId
    }

    # -----------------------------------------------------------------------------
    # Correct the GRanges End Position for INS Calls
    # -----------------------------------------------------------------------------
    # Find all the breakends that belong to an INS call.
    ins_indices <- which(x$type == "INS")
    if (length(ins_indices) > 0) {
      pair_ins_indices <- ins_indices[grepl("_bp2$", names(x)[ins_indices])]
      x[pair_ins_indices] <- IRanges::shift(x[pair_ins_indices], x$insLen[pair_ins_indices])
    }
  }
  return(x)
}

## Check overlaps between a breakpointRanges and a GRanges, 
pairInBed = function(query, subject) {
  
  overlaps = rep(NA, length(query))
  
  ## For each breakend
  processed = c()
  for (i in 1:length(query)) {
    
    partner.idx = which(names(query) == query$partner[i])
    is.translocation = as.character(seqnames(query[i])) != as.character(seqnames(query[partner.idx]))
    
    ## Check if we've already procesed the full breakpoint
    if (i %in% processed) {
      next
    }    
    
    if (is.translocation) {
      
      overlap = any(query[c(i, partner.idx)] %^% subject)      
      
    } else {
      
      bkpt = GRanges(as.character(seqnames(query))[i], IRanges(start(query)[i], end(query)[partner.idx]))  

      # MWL NOTE: the original script check for ANY overlap between the database and the range. 
      #           This was only done for BED files. When BEDPE files are used, 0.8 prop overlap 
      #           is required to consider an overlap valid. 
      #           The code below takes the overlap between SV, and subject DB, 
      #           then for overlapping hits calculates the prop overlap. 
      #           Overlap = TRUE for hits that overlap, and are greater than 80% size overlap.
      hits <- findOverlaps(bkpt, subject)

      if (sum(countSubjectHits(hits)) > 0) {
        hit_overlaps <- pintersect(bkpt[queryHits(hits)], subject[subjectHits(hits)])
        percent_hit_Overlap <- width(hit_overlaps) / width(bkpt[queryHits(hits)])
        if (sum(countSubjectHits(hits[percent_hit_Overlap > 0.8])) > 0) {
          overlap = TRUE
        } else {
          overlap = FALSE
        }
      } else {
        overlap = FALSE
      }
    }
    
    overlaps[c(i, partner.idx)] = overlap
    processed = c(processed, i, partner.idx)
    
  }
  
  return(overlaps)
  
}



## Annotate breakpointRanges object with breakpointRanges or GRanges
annotateDB = function(x, db, name, slop, ignore.strand=F) {
  
  
  ## Use different overlap method depending on whether DB is BED or BEDPE
  if ('partner' %in% colnames(mcols(db))) {
    overlaps = rep(FALSE, length(x))

    # Get unique types present in both x and db
    common_types = intersect(unique(x$type), unique(db$type))

    ### Ensure that type matches between query and db
    ### `findBreakpointOverlaps` struggles with INS calls when sizemargin is used. 
    ### Therefore, a different treatment of INS calls is used. 
    for (svtype in common_types) {
      is_query_type = x$type == svtype
      is_db_type = db$type == svtype

      if (any(is_query_type) && any(is_db_type)) {
      # For INS, do not use sizemargin
      if (svtype == "INS") {
        ov = StructuralVariantAnnotation::findBreakpointOverlaps(
        query = x[is_query_type],
        subject = db[is_db_type],
        maxgap = slop,
        restrictMarginToSizeMultiple = 0.8,
        ignore.strand = ignore.strand
        )
      } else {
        ov = StructuralVariantAnnotation::findBreakpointOverlaps(
        query = x[is_query_type],
        subject = db[is_db_type],
        maxgap = slop,
        sizemargin = 0.8,
        restrictMarginToSizeMultiple = 0.8,
        ignore.strand = ignore.strand
        )
      }
      overlaps[which(is_query_type)[queryHits(ov)]] = TRUE
      }
    }

  } else {
    # overlaps = GenomicRanges::findOverlaps(query=x, subject=db)
    overlaps = pairInBed(query=x, subject=db)
  }
  
  ## Annotate with hits if there are any
  x$db[overlaps] = paste0(dta$db[overlaps],name,',')
  
  return(x)
  
}



## Convert breakpointRanges to BEDPE
vcfToBedpe = function(vcf) {
  
  sqn = as.character(seqnames(vcf))
  strand = as.character(strand(vcf))
  res = c()
  processed = c()
  
  for (i in 1:length(vcf)) {
    bnd = names(vcf)[i]
    partner = vcf$partner[i]
    partner.idx = which(names(vcf) == partner)
    
    ## If we don't have exactly one partner, exclude this variant
    if (length(partner.idx) != 1) {
      warning('Missing partner for breakend ', bnd)
      next
    }
    
    ## Check to see if we've alrady processed this or it's partner
    if (any(c(bnd, partner) %in% processed)) {
      next
    }
    
    ## Aggregate database string
    dbstr = unique(unlist(strsplit(c(vcf$db[i], vcf$db[partner.idx]), ',', fixed=T)))
    dbstr = paste(dbstr, collapse=',')
    dbstr = paste0('known=',dbstr,';')
 
    ## Combine breakends in single line
    res.i = c(sqn[i], start(vcf)[i], end(vcf)[i],                                  ## chr1, start1, end1
              sqn[partner.idx], start(vcf)[partner.idx], end(vcf)[partner.idx],    ## chr2, start2, end 2
              vcf$type[i], '.', strand[i], strand[partner.idx],                    ## type, score, strand1, strand2
              vcf$evidence[i], vcf$tools[i], vcf$SUPP[i], vcf$SUPP_VEC[i], vcf$sampleID[i], dbstr)               ## evidence, tools, TN, info 
    
    ## Add to result, keep track of processed breakends
    res = rbind(res, res.i)
    processed = c(processed, bnd, partner)
  } 

  ## Add colnames and fill in simple event classifications
  colnames(res) = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'type', 'score', 'strand1', 'strand2', 'evidence', 'tools', 'support', 'support_vector', 'sampleID', 'info')
  res = as.data.frame(res, stringsAsFactors=F)

  ## Fix coordinates (have to subtract when starting from a bedpe)
  res$start1 = as.numeric(res$start1) - 1
  res$start2 = as.numeric(res$start2) - 1
  
  
  colnames(res)[1] = paste0('#', colnames(res)[1])
  
  return(res)
  
}



## Collect arguments
option_list = list(
  make_option(c("-b", "--bedpe"),    type='character', help="Input BEDPE"),
  make_option(c("-n", "--db_names"), type='character', help="Comma-delimited list of database names corresponding to the order in --db_files"),
  make_option(c("-f", "--db_files"), type='character', help="Comma-delimited list of database files corresponding to the order in --db_names"),
  make_option(c("-i", "--db_ignore_strand"), type='character', help="Comma-delimited list of database names to ignore strand orientation for when overlapping? Should be present in --db_names"),
  make_option(c("-s", "--slop"),     type='numeric',   help="Padding to use when comparing breakpoints"),
  make_option(c("-o", "--out_file"), type='character', help="Output BEDPE"))
opt = parse_args(OptionParser(option_list=option_list))


## Unpack arguments
opt$db_names = unlist(strsplit(opt$db_names, ',', fixed=T))
opt$db_files = unlist(strsplit(opt$db_files, ',', fixed=T))

if (!is.null(opt$db_ignore_strand)) {
  opt$db_ignore_strand = unlist(strsplit(opt$db_ignore_strand, ',', fixed=T))  
} 



## Sanity-check ignore-strand option
if (!is.null(opt$db_ignore_strand) && !all(opt$db_ignore_strand %in% opt$db_names)) {
  missing = paste(setdiff(opt$db_ignore_strand, opt$db_names), collapse=',')
  stop('Databases present in --db_ignore_strand not present in --db_names: ', missing)
}

## Read bedpe

dta <- tryCatch( 
  {
      readBEDPE(opt$bedpe)
  },
  error = function(e) {
      res <- as.data.frame(matrix(ncol = 16, nrow = 0))
      colnames(res) <- c('#chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'type', 'score', 'strand1', 'strand2', 'evidence', 'tools', 'support', 'support_vector', 'sampleID', 'info')
      colnames(res)[1] <- paste0('#', colnames(res)[1])
      
      write.table(res, opt$out_file, row.names=F, col.names=T, sep='\t', quote=F)
      quit(save = "no", status = 0)
  }
)
# if SV CSV is empty (i.e., no somatic SV), write an empty file as output object.

dta$db = ''

## Annotate with databases
for (i in 1:length(opt$db_names)) {
  
  db.name = opt$db_names[i]
  db.file = opt$db_files[i]
  is = !is.null(opt$db_ignore_strand) && db.name %in% opt$db_ignore_strand

  print(db.name)

  db = readDB(db.file)

  dta = annotateDB(x=dta, 
                   db=db, 
                   name=db.name, 
                   slop=opt$slop,
                   ignore.strand=is)

}

## Convert to bedpe 
res = vcfToBedpe(dta)

## Write result
write.table(res, opt$out_file, row.names=F, col.names=T, sep='\t', quote=F)
