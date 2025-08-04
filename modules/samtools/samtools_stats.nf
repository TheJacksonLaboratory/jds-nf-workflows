process SAMTOOLS_STATS {
    tag "$sampleID"

    cpus 1
    memory 8.GB
    time '06:00:00'

    container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

    publishDir {
        def type = "${params.workflow}" == 'chipseq' ? ( sampleID =~ /INPUT/ ? 'control_samples/' : 'immuno_precip_samples/') : ''
        "${params.pubdir}/${type + sampleID + '/samtools'}"
    }, pattern: "*.flagstat", mode: 'copy'

    publishDir {
        def type = "${params.workflow}" == 'chipseq' ? ( sampleID =~ /INPUT/ ? 'control_samples/' : 'immuno_precip_samples/') : ''
        "${params.pubdir}/${type + sampleID + '/samtools'}"
    }, pattern: "*.idxstats", mode: 'copy'

    publishDir {
        def type = "${params.workflow}" == 'chipseq' ? ( sampleID =~ /INPUT/ ? 'control_samples/' : 'immuno_precip_samples/') : ''
        "${params.pubdir}/${type + sampleID + '/samtools'}"
    }, pattern: "*.stats", mode: 'copy'


    input:
    tuple val(sampleID), file(bam)

    output:
    tuple val(sampleID), file("*.flagstat"), emit: flagstat
    tuple val(sampleID), file("*.idxstats"), emit: idxstat
    tuple val(sampleID), file("*.stats"), emit: stats

    script:

    def prefix = params.workflow == 'chipseq' ? bam[0].baseName : sampleID

    """
    samtools flagstat ${bam[0]} > ${prefix}.flagstat
    samtools idxstats ${bam[0]} > ${prefix}.idxstats
    samtools stats ${bam[0]} > ${prefix}.stats
    """
}
