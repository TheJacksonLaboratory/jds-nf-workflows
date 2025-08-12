process PICARD_MARKDUPLICATES {
    tag "$sampleID"

    cpus 1
    memory 130.GB
    time { bam.size() < 60.GB ? '18:00:00' : '48:00:00' }
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'
    
    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*.txt", mode:'copy'

    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), path("*_dedup.bam"), emit: dedup_bam
    tuple val(sampleID), path("*.txt"), emit: dedup_metrics

    script:
    String my_mem = (task.memory-10.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    picard -Xmx${my_mem}G -Xms${my_mem}G MarkDuplicates \
    I=${bam} \
    O=${bam.baseName}_dedup.bam \
    M=${sampleID}_dup_metrics.txt \
    VALIDATION_STRINGENCY=SILENT \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
    ASSUME_SORT_ORDER="queryname" \
    CLEAR_DT="false" \
    ADD_PG_TAG_TO_READS=false
    """
}
