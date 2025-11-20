process QUILT_TO_QTL2 {

  time {60.min * task.attempt}
  cpus 1
  memory {100.GB * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  
  container 'docker://sjwidmay/bcftools_qtl2:latest'

  publishDir "${params.pubdir}/quilt_vcfs", pattern:"*_merged.vcf.gz", mode:'copy'
  publishDir "${params.pubdir}/quilt_vcfs", pattern:"*_resolution_summary.csv", mode:'copy'

  input:
  tuple val(chr), val(downsample_to_cov), val(start), val(stop), file(sample_genos), file(sample_genos_index), path(covar_file)

  output:
  tuple val(chr), val(downsample_to_cov), path("*_founder_geno.csv"), path("*_sample_geno.csv"), path("*_pmap.csv"), path("*_gmap.csv"), path("covar.csv"), path("pheno.csv"), emit: qtl2files
  tuple val(chr), path("*_resolution_summary.csv"), emit: resolution_summary
  tuple val(chr), path("*_merged.vcf.gz"), emit: imputed_snps

  script:

  """
  # make founder geno table
  bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' ${params.ref_haps_dir}/${params.cross_type}/chr${chr}_phased_snps.vcf.gz | \
        sed 's/[[# 0-9]*]//g' | \
        sed 's/:GT//g' > chr${chr}_fg.txt

  # concatenate chromosome chunks
  bcftools concat ${sample_genos} --output-type z --output chr${chr}_merged.vcf.gz

  # make sample geno table
  bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/HWE\\t%INFO/INFO_SCORE[\\t%GT]\\n' chr${chr}_merged.vcf.gz | \
        sed 's/[[# 0-9]*]//g' | \
        sed 's/:GT//g' > chr${chr}_sg.txt
  
  # quality filter variants and anchor to grid
  Rscript --vanilla ${projectDir}/bin/lcwgs/prepare_qtl2_files.R chr${chr}_fg.txt \
	    chr${chr}_sg.txt \
	    ${covar_file} \
	    ${params.cross_type} \
	    ${params.ref_haps_dir}/${params.cross_type}/chr${chr}_gen_map.txt \
        ${chr} \
        ${params.gridfile}
  """
}
