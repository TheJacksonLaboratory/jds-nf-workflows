process PUBLISH_BAM {
    tag "$sampleID"

    cpus 1
    memory 4.GB
    time '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir "${params.pubdir}/${sampleID + '/bam'}", pattern:"*.bam", mode:'copy', enabled: true
    publishDir "${params.pubdir}/${sampleID + '/bam'}", pattern:"*.bai", mode:'copy', enabled: true

    input:
    tuple val(sampleID), path(bam), path(bai)

    output:
    tuple val(sampleID), path("*.bam"), path("*.bai"), emit: bam_bai

    script:
        """
        mv ${bam} ${sampleID}.bam
        mv ${bai} ${sampleID}.bai
        """
}
