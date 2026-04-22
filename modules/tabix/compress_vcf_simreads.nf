process COMPRESS_INDEX_VCF_SIMREADS {

    cpus = 1
    memory = 6.GB
    time = '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/tabix:1.11--hdfd78af_0'

    input:
    path(vcf)

    output:
    path("*.vcf.gz"), emit: vcf
    path("*.vcf.gz.tbi"), emit: tbi

    script:
    """
    bgzip \
     -c \
     ${vcf} \
     > ${vcf}.gz

    tabix ${vcf}.gz
    """
}
