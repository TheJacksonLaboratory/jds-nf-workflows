process GENOPROBS {

  cpus 4
  memory 360.GB
  time {3.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'docker://sjwidmay/lcgbs_hr:latest'

  input:
  tuple val(project_id), path(cross)

  output:
  tuple val(project_id), file(cross), file("*.rds"), emit: genoprobs

  script:
  
  """
  current_dir=\$(echo pwd)
  hash=\$(\$current_dir | tail -c 9)
  echo \$hash
  Rscript --vanilla ${projectDir}/bin/qtl/calcGenoProbs.R
  mv pr_36state.rds pr_36state_\$hash.rds

  """
}
