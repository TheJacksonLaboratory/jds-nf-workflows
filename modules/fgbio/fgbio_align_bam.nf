process FGBIO_ALIGN_BAM {
    tag "$sampleID"

    cpus 12
    memory 65.GB
    time 48.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'community.wave.seqera.io/library/bwa-mem2_fgbio_samtools:d6fd27126a192efa'

    input:
    tuple val(sampleID), path(unmapped_bam)
    val sort_order // "template-coordinate" (pre-consensus, required by GroupReadsByUmi) or "coordinate" (post-consensus, final)

    output:
    tuple val(sampleID), path("*.mapped.bam"), emit: bam
    tuple val(sampleID), path("*.mapped.bam.bai"), emit: bai, optional: true

    script:
    // samtools fastq -> bwa-mem2 mem -> fgbio ZipperBams, piped: re-aligns the
    // reads while reattaching the UMI tags fgbio ExtractUmis/FastqToBam
    // stored on the unmapped BAM (bypasses Picard's older, separate
    // MergeBamAlignment step -- see notes/milestones.md open item 6).
    def sort_cmd = sort_order == "template-coordinate" ?
        "samtools sort -@ ${task.cpus} --template-coordinate -o ${sampleID}.mapped.bam" :
        "samtools sort -@ ${task.cpus} --write-index -o ${sampleID}.mapped.bam##idx##${sampleID}.mapped.bam.bai"
    """
    samtools fastq -@ ${task.cpus} ${unmapped_bam} \
    | bwa-mem2 mem -t ${task.cpus} -p -K 150000000 -Y ${params.ref_fa_bwamem2_indices} - \
    | fgbio -Xmx${task.memory.giga}g --async-io=true ZipperBams \
        --unmapped ${unmapped_bam} \
        --ref ${params.ref_fa} \
        --output /dev/stdout \
        --tags-to-reverse Consensus \
        --tags-to-revcomp Consensus \
    | ${sort_cmd} -
    """
}
