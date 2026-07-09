process MULTIQC {
    cpus 1
    memory 8.GB
    time '04:00:00'
    
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/multiqc:v1.34.dev0_custom'

    publishDir path: { "${params.pubdir}/multiqc" }, pattern: "*multiqc_report.html", mode:'copy', enabled: !params.save_multiqc_inputs
    publishDir path: { "${params.pubdir}/multiqc" }, pattern: "*_data", mode:'copy', enabled: !params.save_multiqc_inputs
    publishDir path: { "${params.pubdir}/multiqc" }, pattern: "*_plots", mode:'copy', enabled: !params.save_multiqc_inputs
    publishDir path: { "${params.pubdir}/multiqc" }, mode:'copy', enabled: params.save_multiqc_inputs

    input:
    path(multiqc_files)
    path(multiqc_config)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data" , emit: data
    path "*_plots" , optional:true, emit: plots
    path(multiqc_files), emit: input_files

    script:

    def custom_config = multiqc_config ? "--config $multiqc_config" : params.multiqc_config ? "--config ${params.multiqc_config}" : ''
    """
    ln -s ${projectDir}/bin/shared/multiqc/JAX_logo_rgb_transparentback.png .
    ln -s ${projectDir}/bin/shared/multiqc/JAX_logo_white_transparentback.png .
    multiqc . --no-ai ${custom_config}
    """

}

/* NOTE on the configuration file input: 
    This module requires two inputs, The first being a list of files, the second being an optional configuration file. 
    The second channel is either provided with a path to a config, OR an empty list. 
    In 25.10.+ there is the option to accept null for the second channel, but 25.10.+ is not currently accepted (Mar.26)
    
    If provided a config path as input, that will be used.
    If provided an empty list, but a multiqc_config parameter is set, that param will be used.
    If provided an empty list, and no multiqc_config parameter is set, multiqc will run with default settings.
*/