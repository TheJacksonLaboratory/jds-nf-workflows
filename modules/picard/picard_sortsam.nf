process PICARD_SORTSAM {
    tag "$sampleID"

    cpus 1
    memory { sam.size() < 60.GB ? 30.GB : 60.GB }
    time { sam.size() < 60.GB ? '10:00:00' : '24:00:00' }
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

    publishDir "${params.pubdir}/${sampleID + '/bam'}", pattern: "*_sortsam.bam", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), file(sam)
    val(sort_order)

    output:
    tuple val(sampleID), file("*_sortsam.bam"), emit: bam
    tuple val(sampleID), file("*_sortsam.bai"), emit: bai, optional: true

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    /**
     * Determines whether to create an index file during SAM reordering.
     * 
     * Sets the Picard CREATE_INDEX parameter based on the skip_index flag:
     * - If params.skip_index is true: CREATE_INDEX=false (index creation is disabled)
     * - If params.skip_index is false: CREATE_INDEX=true (index creation is enabled)
     * - If the sort order is coordinate and skip_index is not set, index creation is enabled by default since it's generally recommended for coordinate-sorted BAM files. For other sort orders, index creation is disabled by default unless skip_index is explicitly set to false.
     * 
     * Usage: Pass --skip_index to the nextflow command to disable index creation
     */
     
    index_creation = params.skip_index ? 'CREATE_INDEX=false' : (sort_order == 'coordinate' ? 'CREATE_INDEX=true' : '')

    """
    picard -Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp SortSam \
    SO=${sort_order} \
    INPUT=${sam} \
    OUTPUT=${sam.baseName}_sortsam.bam \
    TMP_DIR=`pwd`/tmp \
    VALIDATION_STRINGENCY=SILENT \
    ${index_creation}
    """
}
