process CONCATENATE_SIMREADS {
    tag "$sampleID"

    cpus 1
    memory 15.GB
    time '03:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'ubuntu:20.04'

    publishDir "${params.pubdir}/results/simulated_10x_AllChr", pattern: "*.fastq.gz", mode:'copy', enabled: params.workflow == 'generate_wgs_simreads'
    publishDir "${params.pubdir}/results/simulated_60x_AllChr", pattern: "*.fastq.gz", mode:'copy', enabled: params.workflow == 'generate_wes_simreads'


    input:
    tuple val(sampleID), file(R1), file(R2)

    output:
    tuple val(sampleID), file("*.fastq.gz"), emit: concat_fastq

    script:

    """
    cat $R1 > ${sampleID}_R1${params.extension}
    cat $R2 > ${sampleID}_R2${params.extension}
    """
}
