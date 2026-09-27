process FGBIO_FASTQTOBAM {
    tag "$sampleID"

    cpus 2
    memory 8.GB
    time '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'community.wave.seqera.io/library/bwa-mem2_fgbio_samtools:d6fd27126a192efa'

    input:
    tuple val(sampleID), path(fq_reads)

    output:
    tuple val(sampleID), path("*.unmapped.bam"), emit: bam

    script:
    // Twist UMI adapters: 5bp inline UMI + 2bp skip on both mates, before the
    // template sequence (DOC-001337/DOC-001506). Overridable since a future
    // panel could use a different UMI chemistry.
    """
    fgbio -Xmx${task.memory.giga}g --tmp-dir=. \
    FastqToBam \
    --input ${fq_reads[0]} ${fq_reads[1]} \
    --output ${sampleID}.unmapped.bam \
    --read-structures ${params.umi_read_structure} \
    --sample ${sampleID} \
    --library ${sampleID}
    """
}
