process EDGER_QC {
    tag "$sampleID"

    cpus 6
    memory 85.GB
    time '24:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/biocontainers/mulled-v2-419bd7f10b2b902489ac63bbaafc7db76f8e0ae1:709335c37934db1b481054cbec637c6e5b5971cb-0'

    publishDir "${params.pubdir}/edger", pattern: "hairpin*", mode:'copy'
    publishDir "${params.pubdir}/edger", pattern: "mature*", mode:'copy'


    input:
    path input_files

    output:
    path '*.{txt,pdf,csv}', emit: edger_files


    script:
    """
    ${projectDir}/bin/smrnaseq/edgeR_miRBase.r $input_files
    """

}
