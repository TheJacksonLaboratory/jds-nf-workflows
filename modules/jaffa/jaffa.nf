process JAFFA {

    tag "$sampleID"

    cpus 12
    memory 115.GB
    time 60.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'davidsongroup/jaffa:2.5'

    publishDir path: { "${params.pubdir}/${sampleID + '/fusions'}" }, pattern: "*_jaffa_fusions.csv", mode:'copy'
    publishDir path: { "${params.pubdir}/${sampleID + '/fusions'}" }, pattern: "*_jaffa_fusions.fasta", mode:'copy', enabled: params.keep_intermediate

    input:
        tuple val(sampleID), path(reads)

    output:
        tuple val(sampleID), path("*_jaffa_fusions.csv"), emit: jaffa_fusions
        tuple val(sampleID), path("*_jaffa_fusions.fasta"), emit: jaffa_fasta

    script:
    ext = reads[0].getExtension()

    """

    bpipe run -v \
    -n ${task.cpus} \
    -p fastqInputFormat='*.${ext}' \
    -p refBase=${params.jaffa_ref_dir} \
    -p genome=hg38 \
    -p annotation=genCode22 \
    /JAFFA/JAFFA_direct.groovy \
    ${reads[0]} \
    ${reads[1]}

    mv jaffa_results.csv ${sampleID}_jaffa_fusions.csv
    mv jaffa_results.fasta ${sampleID}_jaffa_fusions.fasta ; 

    """
}
