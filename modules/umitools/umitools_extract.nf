process UMITOOLS_EXTRACT {
    tag "$sampleID"

    cpus 1
    memory 64.GB
    time 10.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container "quay.io/biocontainers/umi_tools:1.1.6--py311haab0aaa_0"

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "${sampleID}.umi.log", mode:'copy'

    input:
    tuple val(sampleID), path(reads)

    output:
        tuple val(sampleID), path("*.umi.log"), emit: log
        tuple val(sampleID), path("*.fastq.gz"), emit: umi_fastq

    script:

    if (params.read_type == 'PE') {
        reads = "-I ${reads[0]} --read2-in=${reads[1]} -S ${sampleID}.umi.R1.fastq.gz --read2-out=${sampleID}.umi.R2.fastq.gz"
    } else {
        reads = "-I ${reads[0]} -S ${sampleID}.umi.R1.fastq.gz"
    }

    pattern_1 = params.umitools_bc_pattern ? "--bc-pattern='${params.umitools_bc_pattern}'" : ''
    pattern_2 = params.umitools_bc_pattern2 ? "--bc-pattern2='${params.umitools_bc_pattern2}'" : ''

    """
    umi_tools extract \
        --umi-separator=${params.umi_separator} \
        -L ${sampleID}.umi.log \
        --extract-method=${params.umitools_extract_method} \
        ${pattern_1} \
        ${pattern_2} \
        ${reads} \
        --temp-dir='./tmp'
    """
}
