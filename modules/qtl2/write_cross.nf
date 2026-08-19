process WRITE_CROSS {

  memory 50.GB
  time 1.hour
  errorStrategy 'retry'
  maxRetries 1

  container 'docker://sjwidmay/lcgbs_hr:latest'

  input:
  tuple val(project_id), path(covar), val(cross_type), path(dedup_samples), path(sampleGenos)
  path(consensusFiles)

  output:
  tuple val(project_id), path("*.rds"), emit: cross

  script:
  """
  current_dir=\$(echo pwd)
  hash=\$(\$current_dir | tail -c 9)
  echo \$hash
  Rscript --vanilla ${moduleDir}/bin/writeControlFile.R
  mv preQC_cross.rds preQC_cross_\${hash}.rds
  """
}
