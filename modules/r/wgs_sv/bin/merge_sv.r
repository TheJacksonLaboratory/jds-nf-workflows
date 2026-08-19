## Merge arbitrary number of VCFs, annotate with simple event type 
libs = c('optparse', 'StructuralVariantAnnotation', 'VariantAnnotation', 'rtracklayer', 'stringr')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T, quietly=T)))
options(width=200, scipen=999)

SUPPORTED_CALLERS = c('manta', 'lumpy', 'delly', 'svaba', 'sniffles', 'pbsv')     ## Update this flag when adding support for new callers

## Callers have different names for the same pieces of evidence,
## For now handle each case separately
getReadSupport = function(vcf, caller, supplementary=FALSE, supported_callers=SUPPORTED_CALLERS) {
  
  ## Don't try to process genotype info if we don't know how
  if (!caller %in% supported_callers) {
    stop('Caller ', caller, ' is not currently supported. Supported callers: ', paste(supported_callers, collapse=','))
  }
    
  sample_id = colnames(geno(vcf)[[1]])

  if (caller == 'manta') {
    ## Common info
    sr = geno(vcf)$SR[, sample_id]
    sr = sapply(sr, `[`, 2)
    pe = geno(vcf)$PR[, sample_id]
    pe = sapply(pe, `[`, 2)
    
    ## Supplementary info
    type = paste0(caller,'_SVTYPE=', info(vcf)$SVTYPE) 
    ft = paste0(caller,'_FT=', geno(vcf)$FT[, sample_id])
    gq = paste0(caller,'_GQ=', geno(vcf)$GQ[, sample_id])
    pl = paste0(caller,'_PL=', geno(vcf)$PL[, sample_id])
    gt = paste0(caller,'_GT=', geno(vcf)$GT[, sample_id])
    homseq = paste0(caller,'_HOMSEQ=', info(vcf)$HOMSEQ)
    homlen = paste0(caller,'_HOMLEN=', info(vcf)$HOMLEN)
    svlen = paste0(caller,'_SVLEN=', info(vcf)$SVLEN)
    
    ### duphold specific.
    dhfc = paste0(caller,'_DHFC=', geno(vcf)$GT[, sample_id])
    dhbfc = paste0(caller,'_DHBFC=', geno(vcf)$GT[, sample_id])
    dhffc = paste0(caller,'_DHFFC=', geno(vcf)$GT[, sample_id])
    dhsp = paste0(caller,'_DHSP=', geno(vcf)$GT[, sample_id])

    supp_string = paste(type, ft, gq, pl, gt, homseq, homlen, svlen, dhfc, dhbfc, dhffc, dhsp, sep=',')

    svlen_string = svlen

  } else if (caller == 'lumpy') {
    
    ## Common info
    sr = info(vcf)$SR
    pe = info(vcf)$PE

    ## Supplementary info
    type = paste0(caller,'_SVTYPE=', info(vcf)$SVTYPE) 
    ro = paste0(caller,'_RO=', geno(vcf)$RO[, sample_id])
    ao = paste0(caller,'_AO=', geno(vcf)$AO[, sample_id])
    dp = paste0(caller,'_DP=', geno(vcf)$DP[, sample_id])
    gt = paste0(caller,'_GT=', geno(vcf)$GT[, sample_id])
    svlen = paste0(caller,'_SVLEN=', info(vcf)$SVLEN)

    ### duphold specific.
    dhfc = paste0(caller,'_DHFC=', geno(vcf)$GT[, sample_id])
    dhbfc = paste0(caller,'_DHBFC=', geno(vcf)$GT[, sample_id])
    dhffc = paste0(caller,'_DHFFC=', geno(vcf)$GT[, sample_id])
    dhsp = paste0(caller,'_DHSP=', geno(vcf)$GT[, sample_id])

    supp_string = paste(type, ro, ao, dp, gt, svlen, dhfc, dhbfc, dhffc, dhsp, sep=',')
    
    svlen_string = svlen

  } else if (caller == 'delly') {
    ## Common info
    sr = info(vcf)$SR
    pe = info(vcf)$PE

    ## Supplementary info
    type = paste0(caller,'_SVTYPE=', info(vcf)$SVTYPE) 
    dr = paste0(caller,'_DR=', geno(vcf)$DR[, sample_id])
    dv = paste0(caller,'_DV=', geno(vcf)$DV[, sample_id])
    rr = paste0(caller,'_RR=', geno(vcf)$RR[, sample_id])
    rv = paste0(caller,'_RV=', geno(vcf)$RV[, sample_id])
    gt = paste0(caller,'_GT=', geno(vcf)$GT[, sample_id])
    homlen = paste0(caller,'_HOMLEN=', info(vcf)$HOMLEN)
    svlen = paste0(caller,'_SVLEN=', info(vcf)$SVLEN)

    ### duphold specific.
    dhfc = paste0(caller,'_DHFC=', geno(vcf)$GT[, sample_id])
    dhbfc = paste0(caller,'_DHBFC=', geno(vcf)$GT[, sample_id])
    dhffc = paste0(caller,'_DHFFC=', geno(vcf)$GT[, sample_id])
    dhsp = paste0(caller,'_DHSP=', geno(vcf)$GT[, sample_id])

    supp_string = paste(type, dr, dv, rr, rv, gt, homlen, svlen, dhfc, dhbfc, dhffc, dhsp, sep=',')

    svlen_string = svlen

  } else if (caller == 'svaba') {
    
    ## Common info
    sr = geno(vcf)$SR[, sample_id]
    pe = geno(vcf)$DR[, sample_id]
    
    ## Supplementary info
    type = paste0(caller,'_SVTYPE=', info(vcf)$SVTYPE) 
    ad = paste0(caller,'_AD=', geno(vcf)$AD[, sample_id])
    dp = paste0(caller,'_DP=', geno(vcf)$DP[, sample_id])
    lo = paste0(caller,'_LO=', geno(vcf)$LO[, sample_id])
    gt = paste0(caller,'_GT=', geno(vcf)$GT[, sample_id])
    homseq = paste0(caller,'_HOMSEQ=', info(vcf)$HOMSEQ) 
    span = paste0(caller,'_SPAN=', info(vcf)$SPAN)

    ### duphold specific.
    dhfc = paste0(caller,'_DHFC=', geno(vcf)$GT[, sample_id])
    dhbfc = paste0(caller,'_DHBFC=', geno(vcf)$GT[, sample_id])
    dhffc = paste0(caller,'_DHFFC=', geno(vcf)$GT[, sample_id])
    dhsp = paste0(caller,'_DHSP=', geno(vcf)$GT[, sample_id])

    supp_string = paste(ad, dp, lo, gt, homseq, span, dhfc, dhbfc, dhffc, dhsp, sep=',')
    
    svlen_string = span
    
  } else if (caller == 'sniffles') {
    
    ## Read support info
    rc = info(vcf)$SUPPORT
    
    # Supplementary Info
    type = paste0(caller,'_SVTYPE=', info(vcf)$SVTYPE)
    svlen = paste0(caller,'_SVLEN=', info(vcf)$SVLEN)
    sv_end = paste0(caller,'_END=', info(vcf)$END)
    gt = paste0(caller,'_GT=', geno(vcf)$GT)
    imprecise = paste0(caller,'_SVTYPE=', info(vcf)$IMPRECISE)
    mosaic = paste0(caller,'_SVTYPE=', info(vcf)$MOSAIC)
  
    supp_string = paste(type, svlen, sv_end, gt, imprecise, mosaic, sep=',')
    
    svlen_string = svlen
    
  } else if (caller == 'pbsv') {
    
    ## Read support info
    rc = c(geno(vcf)$DP)
    
    # Supplementary Info
    type = paste0(caller,'_SVTYPE=', info(vcf)$SVTYPE)
    svlen = paste0(caller,'_SVLEN=', info(vcf)$SVLEN)
    gt = paste0(caller,'_GT=', geno(vcf)$GT)
    sv_end = paste0(caller,'_END=', info(vcf)$END)
    imprecise = paste0(caller,'_SVTYPE=', info(vcf)$IMPRECISE)
    mate = paste0(caller,'_SVTYPE=', info(vcf)$MATEID)
    
    supp_string = paste(type, svlen, sv_end, gt, imprecise, mate, sep=',')
    
    svlen_string = svlen
  }


  # If SVLEN or SPAN is not NA, use it as the event length
  event_length = NA
  if ("svlen" %in% ls() && !all(is.na(svlen))) {
    event_length = svlen
  } else if (exists("span") && !all(is.na(span))) {
    event_length = span
  }

  ## Build output string
  
  # short read callers use split reads and PE reads to assess structural variation, support is different
  if(caller %in% c('manta', 'lumpy', 'delly', 'svaba')){
    # Ensure sr, pe, and svlen_string are always character vectors of the correct length
    sr <- if (is.null(sr) || length(sr) == 0 || all(is.na(sr))) 0 else sr
    pe <- if (is.null(pe) || length(pe) == 0 || all(is.na(pe))) 0 else pe
    svlen_string <- if (is.null(svlen_string) || length(svlen_string) == 0) NA else svlen_string
    
    if (supplementary) {
      res = paste0('[',caller,'_SR=',sr,',', caller,'_PE=', pe, ',', supp_string, ']')
    } else {
      res = paste0('[',caller,'_SR=',sr,',', caller,'_PE=', pe, ',', svlen_string, ']')  
    }
  } else {
    # Ensure rc are always character vectors of the correct length
    rc <- if (is.null(rc) || length(rc) == 0 || all(is.na(rc))) 0 else rc
    svlen_string <- if (is.null(svlen_string) || length(svlen_string) == 0) NA else svlen_string
    
    if (supplementary) {
      res = paste0('[',caller,'_RC=',rc,',', supp_string, ']')
    } else {
      res = paste0('[',caller,'_RC=',rc,',', svlen_string, ']')  
    }
  }
  
  
  return(res)
  
}

sumSupport = function(x) {
  sapply(str_extract_all(x, '(?<=\\=)[0-9]+(?=,|\\])'), function(y) sum(as.numeric(y)))
}

removeRedundantBreakpoints = function(x) {
  # Precompute support and multicaller for all breakends
  read.support <- sumSupport(x$support)
  multicaller.support <- grepl('],[', x$support, fixed = TRUE)
  
  # Build a mapping from breakendPosID to indices
  key.list <- strsplit(x$breakendPosID, ',')
  key.to.idx <- list()
  for (i in seq_along(key.list)) {
    for (k in key.list[[i]]) {
      key.to.idx[[k]] <- c(key.to.idx[[k]], i)
    }
  }
  
  # Find duplicate keys
  key.count <- sapply(key.to.idx, length)
  key.dup <- names(key.count[key.count > 1])
  
  if (length(key.dup) == 0) return(x)
  
  idx.rm <- integer(0)
  
  for (k in key.dup) {
    idx <- key.to.idx[[k]]
    xi <- x[idx]
    xi.read.support <- read.support[idx]
    xi.multicaller <- multicaller.support[idx]
    
    # Keep multi-caller support if present
    idx.multi <- which(xi.multicaller)
    if (length(idx.multi) > 0) {
      idx.rm <- c(idx.rm, idx[-idx.multi])
      next
    }
    
    # Keep highest read support
    idx.max <- which(xi.read.support == max(xi.read.support))
    idx.keep <- idx[idx.max]
    idx.rm.tmp <- setdiff(idx, idx.keep)
    
    # If tie, use SV length or position
    if (length(idx.keep) > 1) {
      xi.keep <- xi[idx.max]
      if (all(!is.na(xi.keep$svLen))) {
        idx.keep.final <- idx.keep[which.max(abs(xi.keep$svLen))]
        idx.rm.tmp <- setdiff(idx.keep, idx.keep.final)
      } else if (all(is.na(xi.keep$svLen)) && length(unique(as.character(seqnames(xi.keep)))) == 1) {
        partners <- x[names(x) %in% xi.keep$partner]
        if (length(partners) > 0) {
          partner.keep <- names(partners)[which.max(start(partners))]
          idx.keep.final <- idx.keep[xi.keep$partner == partner.keep]
          idx.rm.tmp <- setdiff(idx.keep, idx.keep.final)
        }
      }
      # Otherwise, keep all tied SVs
    }
    idx.rm <- c(idx.rm, idx.rm.tmp)
  }
  
  idx.rm <- unique(idx.rm)
  if (length(idx.rm) > 0) {
    x <- x[-idx.rm]
    x <- x[names(x) %in% x$partner]
  }
  return(x)
}


## Compute error between query and subject for a hits object
computeError = function(query, subject, hits) {
  qh <- queryHits(hits)
  sh <- subjectHits(hits)

  # Vectorized local error
  start_q <- start(query)[qh]
  end_q   <- end(query)[qh]
  start_s <- start(subject)[sh]
  end_s   <- end(subject)[sh]

  error_local_min <- pmax(0, pmax(start_q, start_s) - pmin(end_q, end_s))
  error_local_max <- pmax(end_s - start_q, end_q - start_s)

  # For remote error, get partner indices
  query_partner_names   <- as.character(query[qh]$partner)
  subject_partner_names <- as.character(subject[sh]$partner)
  query_partner_idx     <- match(query_partner_names, names(query))
  subject_partner_idx   <- match(subject_partner_names, names(subject))

  # Only compute for valid pairs
  valid <- !is.na(query_partner_idx) & !is.na(subject_partner_idx)
  error_remote_min <- rep(NA_real_, length(qh))
  error_remote_max <- rep(NA_real_, length(qh))
  if (any(valid)) {
    start_qp <- start(query)[query_partner_idx[valid]]
    end_qp   <- end(query)[query_partner_idx[valid]]
    start_sp <- start(subject)[subject_partner_idx[valid]]
    end_sp   <- end(subject)[subject_partner_idx[valid]]

    error_remote_min[valid] <- pmax(0, pmax(start_qp, start_sp) - pmin(end_qp, end_sp))
    error_remote_max[valid] <- pmax(end_sp - start_qp, end_qp - start_sp)
  }

  error <- data.frame(
    local  = error_local_min,
    remote = error_remote_min
  )
  return(error)
}

## Take the union of callsets a and b, both breakpointRanges objects
## If multiple hits found in b for a, choose the closest match, measured
## as the mean distance between breakends
mergeCallsets = function(a, b, slop, margin) {
  
  ## Find overlaps
  cat('Finding overlaps\n')
  overlaps = StructuralVariantAnnotation::findBreakpointOverlaps(query=a, 
                                                                 subject=b, 
                                                                 maxgap=slop, 
                                                                 sizemargin=margin, 
                                                                 restrictMarginToSizeMultiple=0.8)
  
  ## If we have any duplicate query hits, choose hit based on match quality
  if(anyDuplicated(queryHits(overlaps))) {
    
    ## Compute local and remote breakend basepair error on matches
    cat('Computing error for duplicate hits\n')
    error = computeError(query=a, subject=b, hits=overlaps)

    ## Get duplicate hits  
    dup.query.hits = table(queryHits(overlaps))
    dup.query.hits = names(dup.query.hits[dup.query.hits > 1])
    
    ## Determine which hits we're removing
    idx.hits.rm = c()
    for (d in dup.query.hits) {
      
      idx.dup.query.hits = which(queryHits(overlaps) %in% d)
      local.error = error$local[idx.dup.query.hits]
      remote.error = error$remote[idx.dup.query.hits]
      mean.error = rowMeans(cbind(local.error, remote.error))
      
      ## Keep the query hit with the smallest mean error
      idx.hits.rm = c(idx.hits.rm, idx.dup.query.hits[which.max(mean.error)])
      
    }
    overlaps = overlaps[-idx.hits.rm]   
  }
  
  ## For matching SVs, merge caller support
  a$support[queryHits(overlaps)] = paste0(a$support[queryHits(overlaps)],',',b$support[subjectHits(overlaps)])
  a$supplemental[queryHits(overlaps)] = paste0(a$supplemental[queryHits(overlaps)],',',b$supplemental[subjectHits(overlaps)])
  a$breakendPosID[queryHits(overlaps)] = paste0(a$breakendPosID[queryHits(overlaps)],',',b$breakendPosID[subjectHits(overlaps)])
  
  ## Pull in non-matching SVs from b
  res = c(a, b[-subjectHits(overlaps)])

  return(res)
  
}

sup_vector <- function(x) {
  sup_string = ''
  for (caller in opt$callers) {
    if (stringr::str_detect(x, caller)) {
      sup_string <- paste0(sup_string, 1)
    } else {
      sup_string <- paste0(sup_string, 0)
    }
  }
  return(sup_string)
}

## Convert breakpointRanges to BEDPE (optimized for speed)
vcfToBedpe = function(vcf, supplemental=FALSE) {
  n <- length(vcf)
  if (n == 0) {
    # Return empty data.frame with correct columns
    res <- data.frame(
      chr1=character(), start1=integer(), end1=integer(),
      chr2=character(), start2=integer(), end2=integer(),
      type=character(), score=character(),
      strand1=character(), strand2=character(),
      evidence=character(), tools=character(),
      SUPP=integer(), SUPP_VEC=character(),
      stringsAsFactors=FALSE
    )
    colnames(res)[1] <- "#chr1"
    return(res)
  }

  sqn <- as.character(seqnames(vcf))
  strand <- as.character(strand(vcf))
  names_vcf <- names(vcf)
  partner <- vcf$partner
  partner_idx <- match(partner, names_vcf)
  valid <- !is.na(partner_idx) & partner_idx > 0 & partner_idx <= n

  # Only process each pair once: i < partner_idx
  process_mask <- valid & seq_len(n) < partner_idx
  idx <- which(process_mask)
  if (length(idx) == 0) {
    res <- data.frame(
      chr1=character(), start1=integer(), end1=integer(),
      chr2=character(), start2=integer(), end2=integer(),
      type=character(), score=character(),
      strand1=character(), strand2=character(),
      evidence=character(), tools=character(),
      SUPP=integer(), SUPP_VEC=character(),
      stringsAsFactors=FALSE
    )
    colnames(res)[1] <- "#chr1"
    return(res)
  }

  # Pre-allocate vectors
  chr1 <- sqn[idx]
  start1 <- start(vcf)[idx]
  end1 <- end(vcf)[idx]
  chr2 <- sqn[partner_idx[idx]]
  start2 <- start(vcf)[partner_idx[idx]]
  end2 <- end(vcf)[partner_idx[idx]]
  type <- vcf$svtype[idx]
  score <- rep(".", length(idx))
  strand1 <- strand[idx]
  strand2 <- strand[partner_idx[idx]]
  support <- if (supplemental) vcf$supplemental[idx] else vcf$support[idx]

  # Build data.frame
  res <- data.frame(
    chr1=chr1, start1=start1, end1=end1,
    chr2=chr2, start2=start2, end2=end2,
    type=type, score=score,
    strand1=strand1, strand2=strand2,
    evidence=support,
    stringsAsFactors=FALSE
  )

  # Set TRA type for inter-chromosomal events
  res$type[res$chr1 != res$chr2] <- "TRA"

  # Sort by chromosome and coordinates
  chr_levels <- levels(seqnames(vcf))
  res <- res[order(factor(res$chr1, levels=chr_levels), res$start1, res$end1, decreasing=FALSE), , drop=FALSE]

  # Simplify coordinates
  res$end1 <- as.numeric(res$start1) + 1
  res$end2 <- as.numeric(res$start2) + 1

  # Extract tool info from read support column (vectorized)
  tool_matches <- stringr::str_extract_all(res$evidence, "(?<=\\[)[a-z]+(?=_)")
  res$tools <- vapply(tool_matches, function(x) paste(x, collapse=","), character(1))
  res$SUPP <- lengths(tool_matches)
  res$SUPP_VEC <- vapply(res$tools, sup_vector, character(1))

  colnames(res)[1] <- paste0("#", colnames(res)[1])
  return(res)
}



## Collect arguments
option_list = list(
  make_option(c("-v", "--vcf"),                   type='character', help="Comma-delimited list of breakend notation VCFs"),
  make_option(c("-c", "--callers"),               type='character', help="Comma-delimited list of SV caller names corresponding to the order of VCFs given in --vcf"),
  make_option(c("-n", "--sample_name"),           type='character', help="Sample name"),
  make_option(c("-b", "--build"),                 type='character', help="Genome build"),
  make_option(c("-s", "--slop"),                  type='numeric',   help="Padding to use when comparing breakpoints"),
  make_option(c("-m", "--sizemargin"),            type='numeric',   help="Minimum acceptable size fraction of two merged calls"),
  make_option(c("-l", "--min_sv_length"),         type='numeric',   help="Filter SVs shorter than this length"),
  make_option(c("-a", "--allowed_chr"),           type='character', help="Comma-delimited list of chromosomes to keep"),
  make_option(c("-o", "--out_file"),              type='character', help="Output BEDPE"), 
  make_option(c("-p", "--out_file_supplemental"), type='character', help="Output supplemental BEDPE"))
opt = parse_args(OptionParser(option_list=option_list))

# Check that all required options are not NULL
required_opts <- c("vcf", "callers", "sample_name", "build", "slop", "sizemargin", "min_sv_length", "allowed_chr", "out_file", "out_file_supplemental")
missing_opts <- required_opts[sapply(required_opts, function(optname) is.null(opt[[optname]]))]
if (length(missing_opts) > 0) {
  stop("Missing required options: ", paste(missing_opts, collapse = ", "))
}

## Unpack arguments
opt$vcf = unlist(strsplit(opt$vcf, ',', fixed=T))
opt$callers = unlist(strsplit(opt$callers, ',', fixed=T))
opt$allowed_chr = unlist(strsplit(opt$allowed_chr, ',', fixed=T))

cat('Starting....\n')

## Iteratively merge VCFs
res = NULL
for (i in 1:length(opt$vcf)) {
  ## Read VCF
  caller = opt$callers[i]
  cat(paste('Loading calls from', caller, '\n'))
  vcf = VariantAnnotation::readVcf(opt$vcf[i], genome=opt$build)
  # Filter VCF to contain only allowed chromosomes
  vcf = vcf[seqnames(rowRanges(vcf)) %in% opt$allowed_chr, ]
  ## Get read support

  ## If there are no calls in the vcf (edge case), skip that VCF
  if (length(rowRanges(vcf)) == 0) {
    next
  }

  rowRanges(vcf)$support = getReadSupport(vcf=vcf, caller=caller)
  rowRanges(vcf)$supplemental = getReadSupport(vcf=vcf, caller=caller, supplementary=T )

  ## Convert to breakpointRanges object, don't adjust for CIPOS uncertainty (i.e. keep nominalPosition). 
  ## For Manta, infer missing breakpoint is required as the caller does not insert the recip call in the VCF as the other calls do. 
  vcf = StructuralVariantAnnotation::breakpointRanges(vcf, nominalPosition=T, inferMissingBreakends=T)
  
  if(!caller %in% c("pbsv","sniffles")){
    ## Add breakendPosID for later redundancy checks
    vcf$breakendPosID = paste0('[',caller,'=',as.character(seqnames(vcf)),':',start(vcf),':',strand(vcf),']')
  } else {
    ## Add breakendPosID for later redundancy checks
    vcf$breakendPosID = paste0('[',as.character(seqnames(vcf)),':',start(vcf),':*]')
  }

  ## Overlap if this isn't the first callset
  if (i == 1) {
    res = vcf
  } else {
    cat(paste('Merging callset from', caller, 'with prior callsets\n'))
    res = mergeCallsets(a=res, b=vcf, slop=opt$slop, margin=opt$sizemargin)
  }
}

cat('VCF Processing complete\n')

## If res is empty, create an empty data.frame with the expected columns
if (is.null(res) || length(res) == 0) {
  res <- data.frame(
    chr1 = character(),
    start1 = numeric(),
    end1 = numeric(),
    chr2 = character(),
    start2 = numeric(),
    end2 = numeric(),
    type = character(),
    score = character(),
    strand1 = character(),
    strand2 = character(),
    evidence = character(),
    stringsAsFactors = FALSE
  )
  ## Convert to bedpe, apply some filters 
  for (i in c('main','supplemental')) {
    outfile = ifelse(i=='main', opt$out_file, opt$out_file_supplemental)
    write.table(res, outfile, row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)
  }

} else {

  ## Handle breakpoints with duplicate start or end positions
  # length(res)
  cat('Removing redundant breakpoints\n')

  res = removeRedundantBreakpoints(res)

  cat('Final callset contains', length(res), 'breakpoints\n')
  cat('Writing output files\n')

  ## Convert to bedpe, apply some filters 
  for (i in c('main','supplemental')) {
    
    outfile = ifelse(i=='main', opt$out_file, opt$out_file_supplemental)
    
    ## Convert to BEDPE format
    res.i = vcfToBedpe(res, supplemental=i=='supplemental')
    res.i$sampleID = opt$sample_name
    
    ## Filter non-TRA and non-INS variants for minimum length opt$min_sv_length
    sv.lengths = abs(as.numeric(res.i$start2) - as.numeric(res.i$start1))
    res.i = res.i[res.i$type == 'TRA' | res.i$type == 'INS' | sv.lengths >= opt$min_sv_length, ]
    
    ## Filter SVs not occurring in allowed chromosomes (i.e. autosomes and sex chromosomes)
    ## the back ticks are used here in chr1 because of the leading # in the column name. 
    ## the # is to keep the bedpe file compatible with bedtools.
    res.i = res.i[res.i$`#chr1` %in% opt$allowed_chr & res.i$chr2 %in% opt$allowed_chr, ]
    
    ## Write result
    write.table(res.i, outfile, row.names=F, col.names=T, sep='\t', quote=F)
  }
}
cat('Done\n')