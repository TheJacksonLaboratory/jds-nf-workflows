process MOSDEPTH {

  cpus 1
  memory 10.GB
  time '06:00:00'

  container 'mosdepth_v0.3.3'

  publishDir "${params.pubdir}/${sampleID + '/samtools'}", pattern: "*.mosdepth*", mode:'copy'

  input:
  tuple val(sampleID), file(bam), file(index)

  output:
  tuple val(sampleID), file("*.mosdepth*"), emit: mosdepth
  tuple val(sampleID), file("*.regions*"), emit: regions

  script:
  
  """
  mosdepth -n --by 1000 --threads 2 ${sampleID} ${bam}
  """
}