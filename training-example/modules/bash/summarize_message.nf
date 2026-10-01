process SUMMARIZE_MESSAGE {
    tag "${sampleID}"

    input:
    tuple val(sampleID), path(uppercase_file)

    output:
    tuple val(sampleID), path("${sampleID}.summary.txt"), emit: summary_file

    script:
    """
    echo "SampleID: ${sampleID}" > ${sampleID}.summary.txt
    echo "Word count: \$(wc -w < ${uppercase_file})" >> ${sampleID}.summary.txt
    echo "Message:" >> ${sampleID}.summary.txt
    cat ${uppercase_file} >> ${sampleID}.summary.txt
    """
}