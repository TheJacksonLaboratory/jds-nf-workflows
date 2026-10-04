process SNPSIFT_DBNSFP {
    tag "$sampleID"

    cpus 1
    memory  6.GB
    time '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/snpeff_snpsift_5.1:v5.1d'

    // clonal_hema counterpart of snpeff_snpsift/snpsift_dbnsfp.nf (unchanged).
    // Human only -- called from subworkflows/chip_annotation.nf only when
    // params.dbNSFP is set (mouse: no dbNSFP resource; see clonal_hema.config).
    // No publishDir: intermediate only.

    input:
    tuple val(sampleID), file(vcf)
    val(indel_snp)

    output:
    tuple val(sampleID), file("*.vcf"), emit: vcf

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]
    output_suffix = indel_snp == 'INDEL' ? 'INDEL_snpsift_dbNSFPanno.vcf' : 'SNP_snpsift_dbNSFPanno.vcf'
    """
    java -Xmx${my_mem}G -jar /opt/snpEff/SnpSift.jar \
    dbnsfp -v -db ${params.dbNSFP} -noDownload -a \
    -f SIFT_score,SIFT_pred,Polyphen2_HDIV_score,MutationAssessor_score,phyloP100way_vertebrate,1000Gp3_AF,1000Gp3_AFR_AF,1000Gp3_EUR_AF,1000Gp3_AMR_AF,1000Gp3_EAS_AF,ESP6500_AA_AF,ESP6500_EA_AF,aapos,Ensembl_transcriptid,genename \
    ${vcf} > ${sampleID}_${output_suffix}
    """
}
