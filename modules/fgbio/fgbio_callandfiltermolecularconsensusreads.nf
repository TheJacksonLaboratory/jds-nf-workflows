process FGBIO_CALLANDFILTERMOLECULARCONSENSUSREADS {
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
    // Non-duplex counterpart to fgbio_callandfilterduplexconsensusreads.nf,
    // ported from the same reference pipeline (5baseTAPS): collapses each
    // UMI+position family into one consensus read regardless of strand
    // balance -- no both-strand ("duplex") requirement. --min-reads here is
    // a single integer (fgbio's own CallMolecularConsensusReads --help:
    // "--min-reads=1 is appropriate" to retain singletons while still
    // consensus-calling true duplicate families), unlike the duplex tool's
    // 1-3 value syntax -- params.call_min_reads must be set to a single
    // value (e.g. "1") when params.duplex_seq is false.
    def call_mem = (task.memory.giga / 3).intValue()
    def filter_mem = task.memory.giga - call_mem
    """
    fgbio -Xmx${call_mem}g --tmp-dir=. --async-io=true --compression=1 \
    CallMolecularConsensusReads \
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
