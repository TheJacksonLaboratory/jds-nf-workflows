process CONCATENATE_GENOPROBS {

  cpus 8 
  memory {600.GB * task.attempt}
  time {6.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 2 

  container 'docker://sjwidmay/jds_lcwgs_hr:1.0.0'

  publishDir "${params.pubdir}/geno_probs", pattern:"*.rds", mode:'copy'
  
  input:
  tuple val(chrs), val(downsample_to_cov), path(genoprobs), path(crosses)

  output:
  tuple val(downsample_to_cov), path("complete_pmap.rds"), path("complete_genoprobs.rds"), path("complete_alleleprobs.rds"), path("interp_250k_alleleprobs.rds"), path("interp_250k_genoprobs.rds"), path("grid_pmap.rds"), path("interp_250k_kinship_loco.rds"), path("interp_250k_viterbi.rds"), emit: concat_probs

  script:

  """
  Rscript --vanilla ${projectDir}/bin/lcwgs/concatGenoProbs_lcwgs.R ${params.cross_type} ${params.interp_250k_gridfile} ${projectDir}/bin/qtl/interpolate_genoprobs.R ${task.cpus}
  """
}
