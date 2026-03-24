process GENERATE_SIMULATED_RNA_DATA {

    cpus 2
    memory 40.GB
    time 1.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container '/projects/omics_share/meta/containers/simreads.sif'

    publishDir "${params.pubdir}/results", pattern: '*.gz', mode:'copy'
    publishDir "${params.pubdir}/results", pattern: '*true_expression.txt', mode:'copy'


    input:
    path(wrap_fa)

    output:
    path('*.gz')
    path('*true_expression.txt')

    script:
    """
    Rscript ${projectDir}/bin/generate_rnaseq_simreads/generate_simulated_RNA_data.R ${wrap_fa} 
    """

}
