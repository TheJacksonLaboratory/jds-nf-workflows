process SAMTOOLS_FAIDX_TO_BED {
    tag "Fasta"

    cpus 1
    memory 8.GB
    time '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

    output:
        path("*.bed"), emit: bed

    script:
    """
      samtools faidx ${params.ref_fa}
      awk 'BEGIN {FS="\t"}; {print \$1 FS "0" FS \$2}' ${params.ref_fa}.fai > temp_interval_file.bed
    """
}
