process CALC_PBC_METRICS {
    tag "$sampleID"

    cpus 4
    memory 20.GB    
    time '10:00:00' 
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir path: { "${params.pubdir}/${sampleID + '/stats'}" }, pattern: "*.pbc.qc", mode: 'copy'
    
    container 'quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6'

    input:
    tuple val(sampleID), path(tmp_bams)

    output:
    tuple val(sampleID), path("*.pbc.qc")

    script:
    """
    bash "${moduleDir}/bin/bedtools_calc_pbc_metrics.sh" \
        "${tmp_bams[0]}" \
        "${sampleID}"
    """
}
