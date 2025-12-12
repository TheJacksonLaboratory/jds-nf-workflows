process FASTP_LONG {
    tag "$sampleID"

    cpus 6
    memory 85.GB
    time '24:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/biocontainers/fastplong:0.3.0--h224cc79_0'

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "${sampleID}_fastp*.*", mode:'copy'

    input:
        tuple val(sampleID), path(fq)

    output:
        tuple val(sampleID), path("${sampleID}_fastp_report.html"), emit: quality_html
        tuple val(sampleID), path("${sampleID}_fastp.json"), emit: quality_json
        tuple val(sampleID), path("${sampleID}.trimmed.fastq"), emit: trimmed_fastq

    script:

    """
    fastplong -i ${fq} \
        -o ${sampleID}.trimmed.fastq \
        -q ${params.quality_phred} \
        -u ${params.unqualified_perc} \
        -w ${task.cpus} \
        -j ${sampleID}_fastp.json \
        -h ${sampleID}_fastp_report.html \
        -R "${sampleID} fastp report"
    """
}
