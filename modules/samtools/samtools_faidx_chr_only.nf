process SAMTOOLS_FAIDX_CHR_ONLY {
    tag "Fasta"

    cpus 1
    memory 8.GB
    time '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

    input:
        path(fasta)

    output:
        path("*selectedCHR.fa"), emit: chr_fa

    script:
    """
    samtools faidx ${fasta} ${params.chrom_list} > ${fasta.baseName}.selectedCHR.fa
    """

}
