process GS_TO_QTL2 {
  
  cpus 2
  time 1.hour
  memory 50.GB
  errorStrategy {(task.exitStatus == 1) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}.\n Please check ${projectDir}/results/${project_id}/logs.\n Also, please verify that sample names in metadata match those expected in FinalReport file(s).\n\n"; return 'ignore'}.call() : 'finish'}

  container 'docker://sjwidmay/lcgbs_hr:latest'

  input:
  tuple val(project_id), path(finalreport_files), path(covar_file), val(cross_type)

  output:
  tuple val(project_id), path("*geno*.csv"), emit: sampleGenos
  tuple val(project_id), file("*covar.csv"), val(cross_type), emit: qtl2meta
  tuple val(project_id), path("*int.csv"), emit: qtl2ints
  tuple val(project_id), path("*.fst"), emit: qtl2intsfst


  script:
  """
  Rscript --vanilla ${projectDir}/bin/qtl/geneseek2qtl2.R \
	${params.cc_do_allele_codes} \
	${covar_file} \
	${finalreport_files} \
  ${params.max_pct_missing}

  current_dir=\$(echo pwd)
  hash=\$(\$current_dir | tail -c 9)
  mv intensities.fst intensities_\${hash}.fst
  mv covar.csv \${hash}_covar.csv
  mv chrYint.csv chrY_\${hash}_int.csv
  mv chrXint.csv chrX_\${hash}_int.csv
  """
}
