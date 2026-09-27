process FGBIO_GROUPREADSBYUMI {
    tag "$sampleID"

    cpus 4
    memory 16.GB
    time '12:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'community.wave.seqera.io/library/bwa-mem2_fgbio_samtools:d6fd27126a192efa'

    publishDir path: { "${params.pubdir}/${sampleID + '/stats'}" }, pattern: "*.family-sizes.txt", mode:'copy'

    input:
    tuple val(sampleID), path(mapped_bam)

    output:
    tuple val(sampleID), path("*.grouped.bam"), emit: bam
    tuple val(sampleID), path("*.family-sizes.txt"), emit: histogram

    script:
    // strategy=Paired is required for duplex consensus calling (matches
    // DOC-001337 step 7); edits=1 matches the UMI design's minimum Hamming
    // distance of 2 between adapters (DOC-001506), tolerating one sequencing
    // error per UMI without merging genuinely distinct families.
    """
    fgbio -Xmx${task.memory.giga}g --tmp-dir=. \
    GroupReadsByUmi \
    --input ${mapped_bam} \
    --output ${sampleID}.grouped.bam \
    --strategy ${params.groupreadsbyumi_strategy} \
    --edits ${params.groupreadsbyumi_edits} \
    --raw-tag RX \
    --family-size-histogram ${sampleID}.family-sizes.txt
    """
}
