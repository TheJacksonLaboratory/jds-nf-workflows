process GATK_GENOMICSDBIMPORT {
    tag "$sampleID"

    cpus 5
    memory 150.GB
    time '15:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    input:
    tuple val(sampleID), path(gvcf), path(idx), path(interval)

    output:
    tuple val(sampleID), path("${sampleID}_${interval.baseName}_genomicsdb"), emit: genomicsdb

    script:
    // memory needs to be set explicitly
    String my_mem = (task.memory-5.GB).toString()
    my_mem =  my_mem[0..-4]

    inputs = gvcf.collect { "--variant $it" }.join(' ')

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}g -Xms${my_mem}g -Djava.io.tmpdir=`pwd`/tmp" GenomicsDBImport \
    --batch-size 50 \
    --reader-threads ${task.cpus} \
    ${inputs} \
    --intervals ${interval} \
    --genomicsdb-workspace-path ${sampleID}_${interval.baseName}_genomicsdb
    """
}
