process SQUID {

    tag "$sampleID"

    cpus 1
    memory 30.GB
    time 5.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/squid:v1.5.1'

    input:
        tuple val(sampleID), path(bam), path(chimeric_bam)

    output:
        tuple val(sampleID), path("*sv.txt"), emit: squid_fusions

    script:
    """
    squid -b ${bam} -c ${chimeric_bam} -o ${sampleID}.squid.fusions
    """
}
