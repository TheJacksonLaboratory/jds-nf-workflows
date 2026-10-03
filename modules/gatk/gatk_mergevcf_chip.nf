process GATK_MERGEVCF {
    tag "$sampleID"

    cpus 1
    memory 24.GB
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.4.0.0'

    // clonal_hema counterpart of gatk/gatk_mergevcf.nf. Same MergeVcfs command,
    // then bgzip + tabix (clonal_hema's final-output convention: every
    // published VCF is compressed + indexed, e.g. the filtered VCF and the
    // selected_SNP/INDEL VCFs already are). Always published -- this is the
    // final annotated deliverable for clonal_hema, unlike the WES module where
    // this same MergeVcfs step is one of two intermediate merge passes.

    publishDir path: { "${params.pubdir}/${sampleID}" }, pattern: "*_annotated_final.vcf.gz*", mode:'copy'

    input:
    tuple val(sampleID), file(snp_vcf), file(indel_vcf)
    val(suffix)

    output:
    tuple val(sampleID), file("*_${suffix}.vcf.gz"), emit: vcf
    tuple val(sampleID), file("*_${suffix}.vcf.gz.tbi"), emit: idx

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]
    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" MergeVcfs \
    -R ${params.ref_fa} \
    -I ${snp_vcf} \
    -I ${indel_vcf} \
    -O ${sampleID}_${suffix}.vcf.gz
    """
}
