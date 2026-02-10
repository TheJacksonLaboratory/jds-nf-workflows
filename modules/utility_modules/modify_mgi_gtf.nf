process MODIFY_MGI_GTF {
    tag "$gtf"

    cpus 1
    memory 5.GB 
    time '00:20:00'

    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container "quay.io/jaxcompsci/perl:0.1.0"

    publishDir "${params.pubdir}", pattern: "*.gtf", mode:'copy'

    input:
        path(gtf_tmp)
        path(ref_table)
    output:
        path("${gtf_tmp.getBaseName(2)}.gtf"), emit: gtf
        //path("*.gtf"), emit: gtf

    script:
        """
        perl ${projectDir}/bin/rnaseq/modify_mgi_gtf.pl ${gtf_tmp} ${ref_table} 
        """
}

