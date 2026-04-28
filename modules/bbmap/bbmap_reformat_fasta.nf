process REFORMAT_FASTA {

    cpus 2
    memory 4.GB
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bbmap:39.06--h92535d8_0'

    publishDir "${params.pubdir}/results", pattern: "*LINEWRAP.fa", mode:'copy'

    input:
    path(fa_cds)

    output:
    path("*LINEWRAP.fa"), emit: wrap_fasta

    script:
    fa="\$(echo -e ${fa_cds} | sed 's/.gz//g')" 
    suffix="\$(echo -e ${fa} | sed 's/.fa//g')" 
    """
    echo "fa_cds : ${fa_cds}"
    echo "fa : ${fa}"
    echo "suffix : ${suffix}"

    gunzip ${fa_cds}

    reformat.sh in=${fa} out=${suffix}.LINEWRAP.fa fastawrap=80
    """
}
