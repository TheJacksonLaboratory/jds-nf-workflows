process CHIP_GENE_REPORT {

    cpus 1
    memory 8.GB
    time '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    // Same DAG position as MULTIQC: collects every sample's own
    // annotation_table.txt.gz (SNPSIFT_EXTRACTFIELDS, modules/snpeff_snpsift/
    // snpsift_extractfields_chip.nf -- already verified correct, see
    // notes/test_summary_human.txt H15) and combines them into one
    // interactive, sortable/filterable HTML report across the whole run.
    // Reuses the pipeline's own MultiQC container purely for its Python 3
    // runtime (stdlib only, no multiqc itself is invoked) rather than adding
    // a new image for this one step.
    container 'quay.io/jaxcompsci/multiqc:v1.34.dev0_custom'

    publishDir path: { "${params.pubdir}/gene_report" }, pattern: "*.html", mode: 'copy'

    input:
    path(annotation_tables)
    path(report_script)
    path(report_logo)
    val(species)             // 'human' | 'mouse': picks the column set and wording
    val(germline_resource)   // 'yes' | 'no': no = placeholder/empty germline resource (mouse B6 default)
    path(gene_coverage)      // gene_coverage_all_samples.tsv (CLONAL_HEMA_COVERAGE_SUMMARY): Gene Coverage heatmap
    path(on_target)          // on_target_all_samples.tsv: per-sample coverage summary above the heatmap

    output:
    path "*clonalHema_variant_report.html", emit: report

    script:
    // optional run-level prefix, e.g. GTBH26-MengQ-149.clonalHema_variant_report.html (--project_name)
    def report_name = (params.project_name ? "${params.project_name}." : "") + "clonalHema_variant_report.html"
    // consensus mode shown in the per-sample QC table of the report (2026-11-11)
    def consensus_mode = params.duplex_seq_bool ? "duplex, call_min_reads ${params.call_min_reads}" : "single-strand, call_min_reads ${params.call_min_reads}"
    // report_script is a declared `path` input (not a ${moduleDir} string
    // interpolation) specifically so Nextflow's -resume cache hashes its
    // content: a bare moduleDir reference renders to the same literal path
    // string regardless of the script's own content, so an edit to
    // build_gene_report.py would never invalidate a cached, stale output
    // (found by a real -resume rerun producing the pre-edit HTML, 2026-10-02).
    """
    python3 ${report_script} --consensus-mode "${consensus_mode}" --logo ${report_logo} --species ${species} --germline-resource ${germline_resource} --gene-coverage ${gene_coverage} --on-target ${on_target} ${report_name} ${annotation_tables}
    """
}
