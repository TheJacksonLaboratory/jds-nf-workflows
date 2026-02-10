process CREATE_BAMLIST {

  cpus 1
  memory 15.GB
  time '00:30:00'

  container 'ubuntu:20.04'

  input:
  tuple path(bams), val(bam_paths), val(downsample_to_cov)

  output:
  tuple path('bamlist.txt'), val(downsample_to_cov), emit: bam_list

  script:
  """
  printf "%s\\n" ${bam_paths.join(' ')} > bamlist.txt
  """
}