process CONCAT_GENOPROBS {

  cpus 8
  memory {300.GB * task.attempt}
  time {3.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'docker://sjwidmay/lcgbs_hr:latest'
  
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_genoprobs.rds", mode:'copy'
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_alleleprobs.rds", mode:'copy'
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_cross.rds", mode:'copy'
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_maxmarg.rds", mode:'copy'
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_kinship.rds", mode:'copy'
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_genotyping_errors.rds", mode:'copy'
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_bad_markers.rds", mode:'copy'

  input:
  tuple val(project_id), file(crosses), file(genoprobs)

  output:
  tuple val(project_id), file("*_genoprobs.rds"), file("*_alleleprobs.rds"), file("*_cross.rds"), file("*_pmap.rds"), file("*_gmap.rds") ,file("*_maxmarg.rds"), file("*_kinship.rds"), emit: qtl2_files
  tuple val(project_id), file("*_genotyping_errors.rds"), file("*_bad_markers.rds"), emit: qc_files

  script:
  
  """
  echo ${genoprobs} > probs.txt
  echo ${crosses} > crosses.txt
  Rscript --vanilla ${projectDir}/bin/qtl/concatGenoProbs.R
  mv genoprobs.rds ${project_id}_genoprobs.rds
  mv alleleprobs.rds ${project_id}_alleleprobs.rds
  mv cross.rds ${project_id}_cross.rds
  mv pmap.rds ${project_id}_pmap.rds
  mv gmap.rds ${project_id}_gmap.rds
  mv maxmarg.rds ${project_id}_maxmarg.rds
  mv bad_markers.rds ${project_id}_bad_markers.rds
  mv genotyping_errors.rds ${project_id}_genotyping_errors.rds
  mv kinship.rds ${project_id}_kinship.rds
  
  """
}
