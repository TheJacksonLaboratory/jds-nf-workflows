process UPDATE_FILES {

  cpus 8
  memory { params.remove_markers ? 300.GB * task.attempt : 100.GB * task.attempt }
  time {3.hour * task.attempt}
  //errorStrategy 'retry' 
  //maxRetries 1

  container 'docker://sjwidmay/lcgbs_hr:latest'
  
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_updated_genoprobs.rds", mode:'copy'
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_updated_alleleprobs.rds", mode:'copy'
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_updated_cross.rds", mode:'copy'
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_updated_maxmarg.rds", mode:'copy'
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_updated_kinship.rds", mode:'copy'

  input:
  tuple val(project_id), path(covar_file), path(cross_file), path(genoprobs_file), path(alleleprobs_file), path(viterbi_file), path(kinship_file), path(marker_file), val(cross_type)
                                
  output:
  tuple val(project_id), file("*_updated_genoprobs.rds"), file("*_updated_alleleprobs.rds"), file("*_updated_cross.rds"), file("*_updated_maxmarg.rds"), file("*_updated_kinship.rds"), emit: qtl2_files

  script:
  
  """
  Rscript --vanilla ${projectDir}/bin/qtl/updateGenoProbs.R ${covar_file} \
                ${cross_file} \
                ${genoprobs_file} \
                ${alleleprobs_file} \
                ${viterbi_file} \
                ${kinship_file} \
                ${marker_file} \
                ${cross_type} \
                ${params.remove_markers} \
                ${params.correct_ids}

  mv genoprobs.rds ${project_id}_updated_genoprobs.rds
  mv alleleprobs.rds ${project_id}_updated_alleleprobs.rds
  mv cross.rds ${project_id}_updated_cross.rds
  mv maxmarg.rds ${project_id}_updated_maxmarg.rds
  mv kinship.rds ${project_id}_updated_kinship.rds
  
  """
}
