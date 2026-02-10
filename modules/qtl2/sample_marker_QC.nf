process SAMPLE_MARKER_QC {

  cpus 1
  memory {50.GB * task.attempt}
  time {1.hour * task.attempt}
  errorStrategy 'retry'
  maxRetries 1

  container 'docker://sjwidmay/lcgbs_hr:latest'

  input:
  tuple val(project_id), path(crosses), file(excluded_samples), file(intensities)

  output:
  tuple path("QC_1.RData"), val(project_id), emit: qc_data
  tuple path("working_cross.RData"), val(project_id), val(nsamples), emit: genoprobs_cross

  script:
  """
  echo ${crosses} > cross_files.txt
  echo ${excluded_samples} > excluded_samples_files.txt
  echo ${intensities} > intensity_files.txt
  Rscript --vanilla ${projectDir}/bin/scripts/qtl2/sampleQC.R cross_files.txt excluded_samples_files.txt intensity_files.txt
  """
}
