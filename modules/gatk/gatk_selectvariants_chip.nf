process GATK_SELECTVARIANTS {
    tag "$sampleID"

    cpus 1
    memory  6.GB
    time '03:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    // clonal_hema counterpart of gatk_selectvariants.nf (unchanged). Differences:
    // same GATK version as the chip FilterMutectCalls modules (4.4.0.0), bgzipped +
    // tabix-indexed output, and always published: the SNP/INDEL VCFs are final
    // products for clonal_hema, not intermediates (the shared module publishes them
    // only with params.keep_intermediate).
    container 'broadinstitute/gatk:4.4.0.0'

    publishDir path: { "${params.pubdir}/${sampleID}" }, pattern: "*_selected_*.vcf.gz*", mode:'copy'

    input:
    tuple val(sampleID), path(vcf), path(idx)
    val(indel_snp)
    val(suffix)

    output:
    tuple val(sampleID), path("*_${suffix}.vcf.gz"), emit: vcf
    tuple val(sampleID), path("*_${suffix}.vcf.gz.tbi"), emit: idx

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]
    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" SelectVariants \
    -R ${params.ref_fa} \
    -V ${vcf} \
    -select-type ${indel_snp} \
    -O ${sampleID}_${suffix}.vcf.gz
    """
}
