process SPLIT_FILES {
    tag {"split fasta with faidx"}

    cpus 1
    memory 20.GB
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/pyfaidx:0.9.0.4--pyhdfd78af_0'
    
    input:
    path(chr_fa)

    output:
    path("chunks/*.fa")

    script:
    """
    mkdir -p chunks
    faidx --split-files ${chr_fa}
    rm ${chr_fa}
    mv *.fa chunks   
    """
}
