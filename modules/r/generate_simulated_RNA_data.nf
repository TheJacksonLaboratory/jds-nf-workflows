process GENERATE_SIMULATED_RNA_DATA {

    cpus 2
    memory 40.GB
    time 1.hour
    // errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/rna_simreads:v1.0'

    publishDir path: { "${params.pubdir}/results" }, pattern: '*.gz', mode:'copy'
    publishDir path: { "${params.pubdir}/results" }, pattern: '*true_expression.txt', mode:'copy'

    input:
    tuple path(wrap_fa), val(library_size), val(library_strategy)

    output:
    path('*.gz')
    path('*true_expression.txt')

    script:

    quality_reference = params.quality_reference ? params.quality_reference : 'none'

    """
    Rscript ${projectDir}/bin/generate_rnaseq_simreads/generate_simulated_RNA_data.R ${wrap_fa} ${library_size} ${library_strategy} ${params.read_length} ${params.fragment_length_min} ${params.fragment_length_max} ${params.fragment_length_mean} ${params.fragment_length_sd} ${params.simulate_sequencing_error} ${quality_reference}
    """

}
