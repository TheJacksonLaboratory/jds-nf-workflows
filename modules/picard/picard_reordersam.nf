process PICARD_REORDERSAM {
    tag "$sampleID"

    cpus 1
    memory 70.GB
    time '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

    publishDir "${params.pubdir}/${sampleID + '/bam'}", pattern: "*.bam", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), file(bam)
    val(picard_dict)

    output:
    tuple val(sampleID), file("*.bam"), emit: bam
    tuple val(sampleID), file("*.bai"), emit: bai, optional: true

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    /**
     * Determines whether to create an index file during SAM reordering.
     * 
     * Sets the Picard CREATE_INDEX parameter based on the skip_index flag:
     * - If params.skip_index is true: CREATE_INDEX=false (index creation is disabled)
     * - If params.skip_index is false: CREATE_INDEX=true (index creation is enabled)
     * 
     * Usage: Pass --skip_index to the nextflow command to disable index creation
     */
     
    index_creation = params.skip_index ? 'CREATE_INDEX=false' : 'CREATE_INDEX=true'

    """
    picard -Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp ReorderSam \
    INPUT=${bam} \
    OUTPUT=${sampleID}_genome_bam_with_read_group_reorder.bam \
    TMP_DIR=`pwd`/tmp \
    SEQUENCE_DICTIONARY=${picard_dict} \
    ${index_creation}
    """
}
