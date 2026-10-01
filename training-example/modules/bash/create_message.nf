process CREATE_MESSAGE {
    tag "${sampleID}"

    input:
    tuple val(sampleID), val(message)

    output:
    tuple val(sampleID), path("${sampleID}.message.txt"), emit: message_file

    script:
    """
    bash ${moduleDir}/bin/create_message.sh \
        '${message}' \
        '${sampleID}.message.txt'
    """
}