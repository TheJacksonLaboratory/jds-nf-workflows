process QUILT {

  memory {100.GB * task.attempt}
  time {4.hour * task.attempt}
  cpus 1
  errorStrategy 'retry' 
  maxRetries 1

  
  container 'docker://sjwidmay/quilt-nf:latest'

  input:
  tuple path(bamlist), val(downsample_to_cov), val(chr), val(start), val(stop), path(covar_file)

  output:
  tuple val(chr), val(downsample_to_cov), val(start), val(stop), path("quilt.*.vcf.gz"), path("quilt.*.vcf.gz.tbi"), path(covar_file), emit: quilt_vcf


  script:

  """
  Rscript --vanilla ${projectDir}/bin/lcwgs/run_quilt.R ${bamlist} \
      ${chr} \
      ${covar_file} \
      ${params.cross_type} \
      ${params.ref_haps_dir}/${params.cross_type}/chr${chr}.hap.gz \
      ${params.ref_haps_dir}/${params.cross_type}/chr${chr}.samples \
      ${params.ref_haps_dir}/${params.cross_type}/chr${chr}.legend.gz \
      2000 \
      ${start} \
      ${stop}
  """
}
