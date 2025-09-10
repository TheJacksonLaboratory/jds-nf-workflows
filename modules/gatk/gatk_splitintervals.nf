process GATK_SPLITINTERVALS {
    tag "Fasta"
    
    cpus 1
    memory 5.GB
    time '00:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    input:
    path(bed)

    output:
    path("interval-files-folder/*scattered.interval_list"), emit: intervals

    script:
    // memory needs to be set explicitly
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" SplitIntervals \
    -R ${params.ref_fa} \
    -L ${bed} \
    --scatter-count ${params.interval_count} \
    -O interval-files-folder
    """
}
