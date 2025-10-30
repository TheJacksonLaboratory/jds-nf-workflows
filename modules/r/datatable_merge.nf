process TABLE_MERGE {
    tag "$sampleID"

    cpus 6
    memory 85.GB
    time '24:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/biocontainers/r-data.table:1.12.2'


    input:
    path mirtop

    output:
    path "mirna.tsv"   , emit: mirna_tsv


    script:
    """
    ${projectDir}/bin/smrnaseq/collapse_mirtop.r ${mirtop}

    """
}
