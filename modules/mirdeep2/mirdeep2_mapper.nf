def VERSION = '2.0.1'

process MIRDEEP2_MAPPER {
    tag "$sampleID"

    cpus 6
    memory 85.GB
    time '24:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container   'quay.io/biocontainers/mirdeep2:2.0.1.3--hdfd78af_1' 

    publishDir "${params.pubdir}/${sampleID + '/mirdeep'}", pattern: "${sampleID}_*.*", mode:'copy'

    input:
    tuple val(sampleID), path(reads)
    path index

    output:
    tuple val(sampleID), path('*_collapsed.fa'), path('*reads_vs_refdb.arf'), emit: mirdeep2_inputs


    script:
    def index_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
    """
    mapper.pl \\
    $reads \\
        -e \\
        -h \\
        -i \\
        -j \\
        -m \\
        -p $index_base \\
        -s ${sampleID}_collapsed.fa \\
        -t ${sampleID}_reads_vs_refdb.arf \\
        -o 4

    """
}
