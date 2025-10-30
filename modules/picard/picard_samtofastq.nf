process PICARD_SAMTOFASTQ {
    tag "$sampleID"

    cpus = 4
    memory = 15.GB
    time 15.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), file("*.fastq"), emit: fastq

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    mkdir -p tmp
    picard -Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp SamToFastq \
    INPUT=${bam} \
    FASTQ=${bam.baseName}.fastq \
    INTERLEAVE=true \
    NON_PF=true
    """
}
