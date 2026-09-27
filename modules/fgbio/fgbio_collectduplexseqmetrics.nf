process FGBIO_COLLECTDUPLEXSEQMETRICS {
    tag "$sampleID"

    cpus 4
    memory 16.GB
    time '12:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'community.wave.seqera.io/library/bwa-mem2_fgbio_samtools:d6fd27126a192efa'

    publishDir path: { "${params.pubdir}/${sampleID + '/stats'}" }, pattern: "*.txt", mode:'copy'
    publishDir path: { "${params.pubdir}/${sampleID + '/stats'}" }, pattern: "*.pdf", mode:'copy'

    input:
    tuple val(sampleID), path(grouped_bam)

    output:
    tuple val(sampleID), path("*duplex_seq_metrics*.txt"), emit: metrics
    tuple val(sampleID), path("*duplex_seq_metrics*.pdf"), emit: pdf, optional: true

    script:
    // Same input as FGBIO_GROUPREADSBYUMI (the grouped BAM) -- this is a
    // parallel QC branch, not part of the main calling chain. Feeds
    // bin/shared/duplex_metrics.py's fold-reduction/duplex-yield MultiQC
    // summary table.
    """
    fgbio -Xmx${task.memory.giga}g --tmp-dir=. --async-io=true \
    CollectDuplexSeqMetrics \
    --input ${grouped_bam} \
    --output ${sampleID}.duplex_seq_metrics \
    --duplex-umi-counts=true
    """
}
