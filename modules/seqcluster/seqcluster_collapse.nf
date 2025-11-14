process SEQCLUSTER_SEQUENCES {
    tag "$sampleID"

    cpus 6
    memory 85.GB
    time '08:00:00'

    container 'quay.io/biocontainers/seqcluster:1.2.8--pyh5e36f6f_0'

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), path("final/*.fastq.gz"), emit: collapsed

    script:
    """
    seqcluster collapse -f $reads -m 1 --min_size 15 -o collapsed
    gzip collapsed/*_trimmed.fastq
    mkdir final
    mv collapsed/*.fastq.gz final/${sampleID}.fastq.gz
    """

}
