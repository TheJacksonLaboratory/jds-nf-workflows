process TARGET_COVERAGE_METRICS {
    tag "$sampleID"

    cpus 4
    memory 20.GB    
    time '10:00:00' 
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir path: { "${params.pubdir}/${sampleID + '/stats'}" }, pattern: "*coverage_metrics.txt", mode: 'copy'
    
    container 'quay.io/jaxcompsci/bedtools-python3:2.26.0' // note: version difference over other bedtools modules. The 2.23.0 container was failing to parse bed target file. 

    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), path("*coverage_metrics.txt"), emit: qc_metrics

    script:
    """
    bash "${moduleDir}/bin/bedtools_amplicon_metrics.sh" \
        "${params.target_gatk}" \
        "${bam}" \
        "${sampleID}"
    """
}

/*
Calculations Per: https://sfvideo.blob.core.windows.net/sitefinity/docs/default-source/application-note/primerclip-a-tool-for-trimming-primer-sequences-application-note.pdf?sfvrsn=cf83e107_14
*/
