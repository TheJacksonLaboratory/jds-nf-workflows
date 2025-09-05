process BWA_MEM {
    tag "$sampleID"

    cpus 12
    memory 65.GB
    time 12.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bwakit:0.7.17.dev1--hdfd78af_1'

    input:
    tuple val(sampleID), path(fq_reads), val(type)
    val(mapping_index)

    output:
    tuple val(sampleID), path("*.sam"), emit: sam

    script:
    """
    bwa mem -t $task.cpus -K 100000000 -p -v 3 -t 2 -Y ${mapping_index} $fq_reads > ${sampleID}_${type}.sam
    """
}

// Input: "val(index)" refers to an index value for scattered input, not the required BWA mapping index. 
