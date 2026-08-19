process CHAIN_BAD2UNIQ_READS {
    tag "$sampleID"

    cpus 1
    memory 4.GB
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/samtools_with_bc:1.3.1'

    input:
    tuple val(sampleID), path(bad_reads)

    output:
    tuple val(sampleID), path("ReadName_unique"), emit: uniq_reads
    
    when: params.chain != null

    script:
    // Get unique 'bad read names' from bam file using gatk ValidateSamFile out results
    """
    bash "${moduleDir}/bin/samtools_chain_bad2uniq_reads.sh" \
        "${bad_reads}"
    """
}
