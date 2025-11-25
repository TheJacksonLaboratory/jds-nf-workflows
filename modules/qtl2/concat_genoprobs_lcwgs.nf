process CONCATENATE_GENOPROBS {

  cpus 2
  memory {400.GB * task.attempt}
  time {3.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'jds_lcwgs_hr'

  publishDir "${params.pubdir}/geno_probs", pattern:"*.rds", mode:'copy'
  
  input:
  tuple val(chrs), val(downsample_to_cov), file(genoprobs), file(crosses)

  output:
  tuple val(downsample_to_cov), file("complete_pmap.rds"), file("complete_genoprobs.rds"), file("complete_alleleprobs.rds"), file("interp_250k_alleleprobs.rds"), file("interp_250k_genoprobs.rds"), file("grid_pmap.rds"), emit: concat_probs

  script:

  """
  Rscript --vanilla ${projectDir}/bin/lcwgs/concatGenoProbs_lcwgs.R ${params.cross_type} ${params.interp_250k_gridfile} ${projectDir}/bin/qtl/interpolate_genoprobs.R
  """
}
