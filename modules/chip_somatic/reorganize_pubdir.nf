process REORGANIZE_PUBDIR {

    cpus 1
    memory 1.GB
    time '00:10:00'

    container 'quay.io/biocontainers/multiqc:1.28--pyhdfd78af_0'

    // Runs once, last, after every per-sample file has already been
    // published -- see bin/chip_somatic/reorganize_pubdir.py for what it
    // fixes and why (two shared modules, samtools_stats.nf and
    // picard_collecthsmetrics.nf's raw call, publish some chip_somatic
    // output outside each sample's own pubdir folder; not modified here
    // since they're shared with other workflows). Ordering is enforced by
    // taking the MultiQC report as an (unused-content) input -- MultiQC's
    // own input channel only resolves once every per-sample file feeding it
    // has already been produced and published.
    input:
    val(sample_ids)
    path(multiqc_report)

    script:
    """
    ${projectDir}/bin/chip_somatic/reorganize_pubdir.py ${params.pubdir} ${sample_ids.join(' ')}
    """
}
