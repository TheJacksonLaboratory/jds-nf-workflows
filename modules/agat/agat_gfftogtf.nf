

process AGAT_GFFTOGTF {

    cpus 1
    memory 64.GB
    time 1.h

    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container "quay.io/biocontainers/agat:1.5.1--pl5321hdfd78af_0"

    publishDir "${params.pubdir}", pattern: "*.gtf", mode:'copy'

    input:
    path(gff)

    output:
    path('*.gtf'), emit: gtf


    script:
    """
    agat_convert_sp_gff2gtf.pl --gff ${gff} -o ${gff.baseName}.gtf
    """
}

/*
 ------------------------------------------------------------------------------
|   Another GFF Analysis Toolkit (AGAT) - Version: v1.5.1                      |
|   https://github.com/NBISweden/AGAT                                          |
|   National Bioinformatics Infrastructure Sweden (NBIS) - www.nbis.se         |
 ------------------------------------------------------------------------------

At least 1 parameters is mandatory:
Input gff/gtf file (--gff or --gtf).


Usage:
        agat_convert_sp_gff2gtf.pl --gff infile.gff [ -o outfile ]
        agat_convert_sp_gff2gtf -h

MWL NOTE: "All AGAT's scripts with the sp prefix use the AGAT parser, before to perform any supplementary task. So, it is not necessary to run agat_convert_sp_gxf2gxf.pl prior the use of any other sp script."
*/
