process RSEM_SAM_VALIDATOR {
    tag "$sampleID"

    cpus 1
    memory 30.GB
    time 3.h

    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/rsem_bowtie2_star:0.1.0'

    input:
    tuple val(sampleID), path(bam)

    output:
    tuple val(sampleID), env("VALID"), emit: validation

    script:
    """
    rsem-sam-validator ${bam} > ${sampleID}.validation.txt

    if grep -q "The input file is not valid!" ${sampleID}.validation.txt; then
        export VALID=FALSE
    else
        export VALID=TRUE
    fi 
    """
}
