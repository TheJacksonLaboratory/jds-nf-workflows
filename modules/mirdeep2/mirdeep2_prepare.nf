process MIRDEEP2_PIGZ {
    tag "$sampleID"
    
    cpus 2
    memory 12.GB
    time '06:00:00'


    // TODO maybe create a mulled container and uncompress within mirdeep2_mapper?
    container   'quay.io/biocontainers/bioconvert:0.4.3--py_0'

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), path("*.{fastq,fq}"), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pigz -f -d -p $task.cpus $reads

    """

}
