process DUPLEX_METRICS_MQC {
    tag "$sampleID"

    cpus 1
    memory 2.GB
    time '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/multiqc:1.28--pyhdfd78af_0'

    publishDir path: { "${params.pubdir}/${sampleID + '/stats'}" }, pattern: "*_mqc.tsv", mode:'copy'

    input:
    tuple val(sampleID), path(prededup_flagstat), path(postdedup_flagstat), path(raw_hsmetrics), path(duplex_seq_metrics)

    output:
    tuple val(sampleID), path("*_mqc.tsv"), emit: mqc

    script:
    // bin/shared/duplex_metrics.py locates the specific files it needs
    // (*.duplex_yield_metrics.txt, *.family_sizes.txt) by suffix out of the
    // full CollectDuplexSeqMetrics output set -- see that script's own
    // docstring. raw_hsmetrics is a separate CollectHsMetrics run on the
    // pre-consensus BAM (raised COVERAGE_CAP) -- reports the real raw
    // sequencing depth as designed, kept out of MultiQC's own file feed so
    // it doesn't create a second per-sample row (folded into this table's
    // "Raw depth" columns instead). Computes: raw vs. consensus read
    // counts, fold-reduction, duplex yield %, raw depth, and family-size
    // distribution, as two MultiQC custom-content tables.
    """
    ${projectDir}/bin/shared/duplex_metrics.py \
    ${sampleID} \
    ${prededup_flagstat} \
    ${postdedup_flagstat} \
    ${raw_hsmetrics} \
    ${duplex_seq_metrics}
    """
}
