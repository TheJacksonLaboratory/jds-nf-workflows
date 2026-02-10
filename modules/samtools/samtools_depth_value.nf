process SAMTOOLS_DEPTH_VALUE {
  tag "$sampleID"
  
  cpus 1
  memory 10.GB
  time '3:00:00'

  container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'

  input:
  tuple val(sampleID), file(bam), file(bai)

  output:
  tuple val(sampleID), file(bam), emit: bam_out
  tuple val(sampleID), file("*_GW_coverage.txt"), emit: depth_coef

  script:

  """
  samtools depth ${bam} -a | awk 'BEGIN{sum=0} {sum += \$3; n++} END{print sum/n}' > ${sampleID}_GW_coverage.txt
  """
}