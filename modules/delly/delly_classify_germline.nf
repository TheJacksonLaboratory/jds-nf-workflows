process DELLY_CLASSIFY {
    tag "$sampleID"
    
    cpus = 1
    memory 80.GB
    time '12:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/biocontainers/delly:1.1.6--h6b1aa3f_2'
    
    publishDir "${params.pubdir}/${sampleID + '/callers'}", pattern: "*.bcf", mode: 'copy'

    input:
    tuple val(sampleID), path(bcf), path(csi)

    output:
    tuple val(sampleID), path("*.bcf"), path("*.csi"), emit: bcf_csi

    script:
    """
    delly classify -f germline -o ${sampleID}_delly_cnv_classified.bcf ${bcf}
    """
}
