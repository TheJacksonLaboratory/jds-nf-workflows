def VERSION = '2.0.1'

process MIRDEEP2_RUN {
    tag "$sampleID"

    cpus 6
    memory 85.GB
    time '24:00:00'
    errorStrategy 'ignore'

    container 'quay.io/biocontainers/mirdeep2:2.0.1.3--hdfd78af_1'

    publishDir "${params.pubdir}/${sampleID + '/mirdeep2'}", pattern: "result*.{bed,csv,html}", mode:'copy'

    input:
    path(fasta)
    tuple val(sampleID), path(reads), path(arf)
    path(hairpin)
    path(mature)

    output:
    path 'result*.{bed,csv,html}', emit: result


    script:
    """
    miRDeep2.pl  \\
        $reads   \\
        $fasta   \\
        $arf     \\
        $mature  \\
        none     \\
        $hairpin \\
        -d       \\
        -z _${reads.simpleName}

    """
}

