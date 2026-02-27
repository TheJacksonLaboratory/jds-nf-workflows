process MODIFY_MGI_GTF {
    tag "$gtf"

    cpus 1
    memory 5.GB 
    time '00:20:00'

    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container "quay.io/jaxcompsci/perl:0.1.0"

    publishDir "${params.pubdir}", pattern: "*.gtf", mode:'copy'

    input:
        path(gtf)
        path(ref_table)
    output:
        path("${gtf.getBaseName(2)}.biotypeAdded.gtf"), emit: gtf

    // note: getBaseName(2) is used to remove both suffixes ".converted.gtf" from the output of the previous step.

    script:
        """
        perl ${projectDir}/bin/generate_rnaseq_index/modify_mgi_gtf.pl ${gtf} ${ref_table} ${gtf.getBaseName(2)}.biotypeAdded.gtf
        """
}
