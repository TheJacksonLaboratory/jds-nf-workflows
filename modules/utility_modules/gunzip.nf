process GUNZIP {

    cpus 1  
    memory 5.GB
    time 2.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container "quay.io/jaxcompsci/py3_perl_pylibs:v2"

    input:
    path(file)

    output:
    path("*"), emit: gunzip_file
    
    script:
    """
    gunzip -c ${file} > ${file.baseName}
    """
}
