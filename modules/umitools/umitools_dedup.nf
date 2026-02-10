process UMITOOLS_DEDUP {
    
    tag "$sampleID"

    cpus 1
    memory 64.GB
    time 10.h

    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container "quay.io/biocontainers/umi_tools:1.1.6--py311haab0aaa_0"

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "${sampleID}.dedup_stats*", mode:'copy'
    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "${sampleID}.dedup_log.txt", mode:'copy'

    input:
    tuple val(sampleID), path(bam), path(bai)

    output:
        tuple val(sampleID), path("*.dedup.bam"), emit: bam
        tuple val(sampleID), path("*.dedup_log.txt"), emit: log
        tuple val(sampleID), path("*.dedup_stats*"), emit: stats

    script:
    paired_option = params.read_type == 'PE' ? '--paired' : ''

    """
    umi_tools dedup \
        -I ${bam} \
        -S ${bam.baseName}.dedup.bam \
        --umi-separator=${params.umi_separator} \
        --multimapping-detection-method=NH \
        --output-stats=${sampleID}.dedup_stats \
        --log=${sampleID}.dedup_log.txt \
        ${paired_option}
    """
}
