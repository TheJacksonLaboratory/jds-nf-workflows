process MAKE_VCF_LIST {
    tag "$sampleID"
    
    cpus 1
    time '00:05:00'

    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'ubuntu:20.04'

    input:
    tuple val(sampleID), val(chroms)
    val(order)

    output:
    tuple val(sampleID), file("*.list"), emit: list

    script:
    // Puts Individual Chromosome Files In Order and Then Into List for MergeVCFs
    // convert paths to strings
    def string_list = chroms.collect { it -> it.toString() }
    // find matches and put in final list
    def sorted = ''
    order.each { chr ->
        sorted += (string_list.find { it -> it.contains('_' + chr + '.vcf') }) + "\n"
    }

    """
    echo "$sorted" > ${sampleID}.list
    """
}
