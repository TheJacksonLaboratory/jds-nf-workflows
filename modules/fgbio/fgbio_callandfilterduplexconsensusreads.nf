process FGBIO_CALLANDFILTERDUPLEXCONSENSUSREADS {
    tag "$sampleID"

    cpus 8
    memory 24.GB
    time 24.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'community.wave.seqera.io/library/bwa-mem2_fgbio_samtools:d6fd27126a192efa'

    input:
    tuple val(sampleID), path(grouped_bam)

    output:
    tuple val(sampleID), path("*.cons.filtered.bam"), emit: bam

    script:
    // CallDuplexConsensusReads piped directly into FilterConsensusReads,
    // matching DOC-001337 step 8 exactly. params.call_min_reads defaults to
    // "2 1 1" (Twist's recommendation for error-correction/low-frequency
    // variant calling, DOC-001506) but is left overridable at run time --
    // a 1-1-0 filter is Twist's separate recommendation for pure duplicate
    // *counting* use cases (e.g. RNA-seq), not this one.
    def call_mem = (task.memory.giga / 3).intValue()
    def filter_mem = task.memory.giga - call_mem
    """
    fgbio -Xmx${call_mem}g --tmp-dir=. --async-io=true --compression=1 \
    CallDuplexConsensusReads \
    --input ${grouped_bam} \
    --output /dev/stdout \
    --min-reads ${params.call_min_reads} \
    --min-input-base-quality ${params.call_min_baseq} \
    --threads ${task.cpus} \
    | fgbio -Xmx${filter_mem}g --tmp-dir=. --compression=0 \
    FilterConsensusReads \
    --input /dev/stdin \
    --output ${sampleID}.cons.filtered.bam \
    --ref ${params.ref_fa} \
    --min-reads ${params.call_min_reads} \
    --min-base-quality ${params.call_min_baseq} \
    --max-base-error-rate ${params.filter_max_base_error_rate}
    """
}
