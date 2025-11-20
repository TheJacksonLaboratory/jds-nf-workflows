process QTL2_GENOPROBS {

  cpus 2
  memory {200.GB * task.attempt}
  time {12.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'jds_lcwgs_hr'

  publishDir "${params.pubdir}/geno_probs", pattern:"*.RData", mode:'copy'
  publishDir "${params.pubdir}/geno_probs", pattern:"covar.csv", mode:'copy'
  
  input:
  tuple val(chr), val(downsample_to_cov), path(founder_geno), path(sample_genos), path(pmap), path(gmap), path(covar), path(pheno)

  output:
  tuple val(chr), val(downsample_to_cov), path("*36_state_probs.RData"), path("*_cross.RData"), emit: geno_probs_out

  script:

  """
  Rscript --vanilla ${projectDir}/bin/lcwgs/genoprobs.R ${chr} \
	${sample_genos} \
	${founder_geno} \
	${pmap} \
	${gmap} \
	${covar} \
	${params.cross_type} \
    ${projectDir}/bin/qtl/smooth_genoprobs.R \
    ${params.smooth_window}
  """
}
