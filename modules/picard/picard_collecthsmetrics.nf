process PICARD_COLLECTHSMETRICS {
    tag "$sampleID"

    cpus 1
    memory  6.GB
    time '08:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

    publishDir path: { "${params.pubdir}/${sampleID + '/stats'}" }, pattern: "*.txt", mode:'copy'

    input:
    tuple val(sampleID), file(bam), file(bai)

    output:
    tuple val(sampleID), file("*Metrics.txt"), emit: hsmetrics

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]
    // No-op unless the caller's config defines this param. Picard's own
    // default (200) caps the per-base histogram used for MEDIAN_TARGET_COVERAGE
    // and the PCT_TARGET_BASES_*X columns -- fine at WGS/WES depth, silently
    // truncates those two metrics at panel depths above ~200x.
    coverage_cap = params.containsKey('coverage_cap') ?
        "COVERAGE_CAP=${params.coverage_cap}" : ''

    """
    picard -Xmx${my_mem}G CollectHsMetrics \
    INPUT=${bam} \
    OUTPUT=${sampleID}_CoverageMetrics.txt \
    BAIT_INTERVALS=${params.bait_picard} \
    TARGET_INTERVALS=${params.target_picard} \
    REFERENCE_SEQUENCE=${params.ref_fa} \
    ${coverage_cap} \
    VALIDATION_STRINGENCY=SILENT
    """
}
