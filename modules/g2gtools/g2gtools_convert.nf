process G2GTOOLS_CONVERT {
    tag "$strain"

    cpus 1
    memory 5.GB
    time '02:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/g2gtools:74926ad'

    publishDir path: { "${params.pubdir}/g2gtools" }, pattern: "*.*", mode:'copy', enabled: { params.append_chromosomes == 'false' || !params.append_chromosomes }

    input:
    tuple val(strain), path(vci), path(tbi)
    path(input_file)
    val(format)
    val(reverse)

    output:
    tuple val(strain), path("*.${format}"), emit: coverted_file
    tuple val(strain), path("*.${format}.unmapped"), emit: unmapped_file

    script:

    def debug_run = params.debug ? '--debug' : ''
    def run_reverse = reverse ? '--reverse' : ''

    """
    g2gtools convert ${debug_run} -i ${input_file} -c ${vci} --file-format ${format} ${run_reverse} -o ${strain}.${params.genome_version}.${format}
    """

    stub:
    """
    touch ${strain}.${params.genome_version}.${format}
    """

}

/*
    ----- tool tip ----------
    Usage: g2gtools convert [OPTIONS]

    Convert coordinates of BAM|SAM|GTF|GFF|BED files

    ╭─ Options ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╮
    │ *  --input-file   -i      FILE                   Input file to convert to new coordinates [default: None] [required]                                                                                                                                                                  │
    │ *  --vci-file     -c      FILE                   VCI file [default: None] [required]                                                                                                                                                                                                  │
    │ *  --file-format  -f      [BAM|SAM|GFF|GTF|BED]  Input file format [default: None] [required]                                                                                                                                                                                         │
    │    --output       -o      FILE                   Name of output file [default: None]                                                                                                                                                                                                  │
    │    --reverse      -r                             Reverse the direction of the conversion                                                                                                                                                                                              │
    │    --verbose      -v      INTEGER                specify multiple times for more verbose output [default: 0]                                                                                                                                                                          │
    │    --help                                        Show this message and exit.                                                                                                                                                                                                          │
    ╰───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯

    ----------------------

    Note regarding 'append_chromosomes': Sanger in does not include 'Y' or 'MT' in the VCF file (Y for biological reasons, MT for seq depth reasons). 
                                         Having the option to map to Y and MT should be available. This extra block of code appends grep 
                                         matched strings: "<STRING>'\t'" bottom to the GTF file.
*/
