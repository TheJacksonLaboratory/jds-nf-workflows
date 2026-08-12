process SAMTOOLS_REHEADER_RGSM {
    tag "$sampleID"

    cpus 1
    memory 8.GB
    time '02:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

    publishDir "${params.pubdir}/${sampleID + '/bam'}", pattern:"*.bam", mode:'copy', enabled: params.containsKey('merge_inds') && params.merge_inds

    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), path("*.bam"), emit: bam

    script:
    """
    samtools view -H ${bam} \
      | awk -v sm='${sampleID}' 'BEGIN{FS=OFS="\\t"} /^@RG/ {for(i=1;i<=NF;i++){if(\$i ~ /^SM:/){\$i="SM:"sm}}} {print}' \
      > ${sampleID}.reheader.sam

    samtools reheader ${sampleID}.reheader.sam ${bam} > ${sampleID}_rgsm_normalized.bam
    """
}
