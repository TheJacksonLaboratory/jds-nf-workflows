

process GFFREAD_GFF3TOGTF {

    cpus 1
    memory 64.GB
    time 1.h

    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container "quay.io/biocontainers/gffread:0.12.1--h2e03b76_1"

    //publishDir "${params.pubdir}", pattern: "*.gtf", mode:'copy'

    input:
    path(gff3)

    output:
    path('*tmp.gtf'), emit: gtf


    script:
    """
    gffread ${gff3} -T -o ${gff3.baseName}.tmp.gtf
    
    """
}

/*
 --------------------------------------------------------------------------------------------
|   The GFFREAD Utility (GFFREAD) - Version: v0.12.1                                         |
|   https://github.com/gpertea/gffread                                                       |
|   Johns Hopkins University Center for Computational Biology (CCB) - https://ccb.jhu.edu/   |
 --------------------------------------------------------------------------------------------

At least 1 parameters is mandatory:
Input gff3/gtf file (--gff3 or --gtf).


Usage:
        gffread infile.gff3 -T [ -o outfile ]
        gffread -h

*/
