process UPPERCASE_MESSAGE {
    tag "${sampleID}"

    input:
    tuple val(sampleID), path(message_file)

    output:
    tuple val(sampleID), path("${sampleID}.uppercase.txt"), emit: uppercase_file

    script:
    """
    tr '[:lower:]' '[:upper:]' < ${message_file} > ${sampleID}.uppercase.txt
    """
}