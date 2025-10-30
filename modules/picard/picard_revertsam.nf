process PICARD_REVERTSAM {
    tag "$sampleID"

    cpus = 4
    memory = 15.GB
    time 15.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), file("*.bam"), emit: bam

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    mkdir -p tmp
    picard -Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp RevertSam \
    INPUT=${bam} \
    OUTPUT_BY_READGROUP=false \
    OUTPUT=${bam.baseName}.reverted.bam \
    VALIDATION_STRINGENCY=LENIENT \
    ATTRIBUTE_TO_CLEAR=FT \
    ATTRIBUTE_TO_CLEAR=CO \
    SORT_ORDER=queryname \
    RESTORE_ORIGINAL_QUALITIES=false
    """
}
