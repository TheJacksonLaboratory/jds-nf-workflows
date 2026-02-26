process CONCAT_INTENSITIES {

  cpus 1
  memory 40.GB
  time 1.hour

  container 'docker://sjwidmay/lcgbs_hr:latest'
  
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_all_intensities.fst", mode:'copy'
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*_intensities.csv", mode:'copy'
  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"dedup_samples.csv", mode:'copy'

  input:
  tuple val(project_id), file(sex_chromosome_intensities), file(intensity_fsts)

  output:
  tuple val(project_id), file("*_chrX_intensities.csv"), file("*_chrY_intensities.csv"), file("*_all_intensities.fst"), emit: concat_intensities
  tuple val(project_id), file("dedup_samples.csv"), emit: dedup_samples

  script:
  
  """
  Rscript --vanilla ${projectDir}/bin/qtl/concatIntensities.R
  mv chrX_intensities.csv ${project_id}_chrX_intensities.csv
  mv chrY_intensities.csv ${project_id}_chrY_intensities.csv
  mv all_intensities.fst ${project_id}_all_intensities.fst
  """
}
