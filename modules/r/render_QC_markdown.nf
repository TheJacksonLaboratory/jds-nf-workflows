process QC_REPORT {

  cpus 1
  memory 50.GB
  time 2.hour

  container 'docker://sjwidmay/haplotype_reconstruction_qtl_nf:qc_markdown'

  publishDir "${params.pubdir}/projects/${project_id}/results", pattern:"*", mode:'copy'

  input:
  tuple val(project_id), path(genoprobs), path(alleleprobs), path(cross), path(viterbi), path(kinship), path(genotyping_errors), path(bad_markers), path(x_intensities), path(y_intensities), path(all_marker_intensities)

  output:
  tuple path("sample_QC.csv"), path("QC_markdown.html"), path("QC_markdown.Rmd"), emit: qc_markdown

  script:

  """
  ls ${projectDir}/bin/qtl/QC_template.Rmd
  cat ${projectDir}/bin/qtl/QC_template.Rmd > QC_working.Rmd
  Rscript --vanilla ${projectDir}/bin/qtl/render_markdown.R QC_working.Rmd
  mv QC_working.html QC_markdown.html
  mv QC_working.Rmd QC_markdown.Rmd
  """
}
