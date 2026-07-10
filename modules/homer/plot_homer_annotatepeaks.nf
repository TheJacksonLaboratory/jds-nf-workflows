process PLOT_HOMER_ANNOTATEPEAKS {

    cpus 2
    memory 10.GB
    time '10:00:00'

    container 'quay.io/biocontainers/mulled-v2-ad9dd5f398966bf899ae05f8e7c54d0fb10cdfa7:05678da05b8e5a7a5130e90a9f9a6c585b965afa-0'

    publishDir path: { "${params.pubdir}/immuno_precip_samples/cross_sample_plots" }, pattern: "*.{pdf,txt}", mode: 'copy'

    input:
    path(annos)
    path(mqc_header)
    val(suffix)

    output:
    path('*.txt'), emit: txt
    path('*.pdf'), emit: pdf
    path('*.tsv'), emit: tsv

    when:
    params.macs_gsize && !params.skip_peak_annotation && !params.skip_peak_qc

    // This script was bundled within the nf-core/chipseq/bin/ directory
    script:
    def prefix = "macs_annotatepeaks"
    def joinedAnnos = annos.join(',')
    def sampleNames = annos.collect { it.name.replace(suffix, '') }.join(',')
    """
    ${moduleDir}/bin/plot_homer_annotatepeaks.r \\
        -i ${joinedAnnos} \\
        -s ${sampleNames} \\
        -p $prefix \\
        -o ./

    find ./ -type f -name "*summary.txt" -exec cat {} \\; | cat $mqc_header - > ${prefix}.summary_mqc.tsv

    """
}
